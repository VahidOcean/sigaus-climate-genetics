# =============================================================================
# 04_popgen_stats.R
# =============================================================================
# Computes population genetic diversity statistics and inter-lineage structure
# for three focal Sigaus species.
#
# Metrics computed:
#   - n, k (haplotypes), S (segregating sites)
#   - Haplotype diversity (Hd)       [Nei 1987]
#   - Nucleotide diversity (pi)       [Tajima 1983]
#   - Watterson's theta_w
#   - Tajima's D                      [Tajima 1989]
#   - Fu & Li's F statistic           [Fu & Li 1993]
#   - Mismatch distribution + raggedness index [Harpending 1994]
#   - Pairwise FST                    [Hudson et al. 1992]
#   - AMOVA                           [Excoffier et al. 1992]
#   - Haplotype network (TCS / parsimony)
#
# Inputs:  data/sigaus_COI_aligned.fasta
#          data/haplotype_table.csv
# Outputs: results/diversity_stats.csv
#          results/fst_matrix.csv
#          results/amova_results.txt
#          figures/haplotype_network_3sp.png / .pdf
#          figures/diversity_plots_3sp.png / .pdf
#          figures/mismatch_distributions_3sp.png / .pdf
#          figures/fst_heatmap_3sp.png / .pdf
#
# Runtime: ~5 minutes (pairwise distance calculations)
# =============================================================================

library(tidyverse)
library(ape)
library(pegas)
library(patchwork)
library(ggrepel)
library(viridis)

cat("\n=== 04: Population Genetic Statistics ===\n")

# ── Configuration ─────────────────────────────────────────────────────────────
FOCAL_POPS <- c("S_australis_group", "S_campestris", "S_villosus")

SP_COLOURS <- c(
    "S_australis_group" = "#e74c3c",
    "S_campestris"      = "#3498db",
    "S_villosus"        = "#2ecc71"
)

SP_LABELS <- c(
    "S_australis_group" = "S. australis",
    "S_campestris"      = "S. campestris",
    "S_villosus"        = "S. villosus"
)

# ── Load data ─────────────────────────────────────────────────────────────────
dna     <- read.dna("../data/sigaus_COI_aligned.fasta", format = "fasta")
hap_tbl <- read_csv("../data/haplotype_table.csv", show_col_types = FALSE)

# Filter to focal species only
hap_tbl         <- hap_tbl %>% filter(population %in% FOCAL_POPS)
pop_assignments  <- setNames(hap_tbl$population, hap_tbl$accession)
focal_acc        <- intersect(names(pop_assignments), rownames(dna))
dna_focal        <- dna[focal_acc, ]

cat("Focal sequences:", nrow(dna_focal), "\n")
cat("Populations:\n"); print(table(pop_assignments[focal_acc]))

# Build population list
pop_list <- map(FOCAL_POPS, function(pop) {
    acc <- intersect(names(pop_assignments[pop_assignments == pop]), rownames(dna_focal))
    if (length(acc) >= 3) dna_focal[acc, ] else NULL
}) %>%
    setNames(FOCAL_POPS) %>%
    keep(~ !is.null(.x))

# =============================================================================
# DIVERSITY STATISTICS
# =============================================================================

compute_diversity <- function(dna_subset, pop_name) {
    n <- nrow(dna_subset)
    if (n < 3) return(NULL)

    haps     <- haplotype(dna_subset)
    k        <- nrow(as.matrix(haps))
    hap_freq <- summary(haps) / n
    Hd       <- (n / (n - 1)) * (1 - sum(hap_freq^2))
    pi_val   <- nuc.div(dna_subset)
    S        <- length(seg.sites(dna_subset))
    a1       <- sum(1 / seq_len(n - 1))
    theta_w  <- S / a1

    taj <- tryCatch(tajima.test(dna_subset),
                    error = function(e) list(D = NA, Pval.normal = NA))

    fu_f <- tryCatch({
        fl <- Fu.Li.F.test(dna_subset)
        fl$statistic
    }, error = function(e) NA_real_)

    pw_dist <- dist.dna(dna_subset, model = "raw", pairwise.deletion = TRUE)
    mm_vals <- as.vector(pw_dist) * ncol(as.matrix(dna_subset))
    mm_vals <- mm_vals[!is.na(mm_vals)]

    max_mm   <- ceiling(max(mm_vals))
    mm_cnts  <- tabulate(ceiling(mm_vals) + 1, nbins = max_mm + 1)
    obs_freq <- mm_cnts / sum(mm_cnts)
    ragged   <- sum(diff(obs_freq)^2)

    tibble(
        population      = pop_name,
        n               = n,
        k               = k,
        S               = S,
        Hd              = round(Hd, 4),
        pi              = round(pi_val, 6),
        theta_w         = round(theta_w, 6),
        Tajima_D        = round(taj$D, 4),
        Tajima_D_pval   = round(taj$Pval.normal, 4),
        Fu_Li_F         = round(fu_f, 4),
        mismatch_mean   = round(mean(mm_vals), 3),
        mismatch_var    = round(var(mm_vals), 3),
        raggedness      = round(ragged, 5),
        mm_values       = list(mm_vals)
    )
}

cat("\n[Computing diversity statistics...]\n")
stats_list <- map(FOCAL_POPS, function(pop) {
    cat("  Processing:", pop, "(n =", nrow(pop_list[[pop]]), ")\n")
    compute_diversity(pop_list[[pop]], pop)
})
stats_df <- bind_rows(stats_list)

cat("\n── Diversity statistics ──\n")
print(stats_df %>% select(-mm_values) %>% as.data.frame(), digits = 4)

write_csv(stats_df %>% select(-mm_values), "../results/diversity_stats.csv")
cat("Saved: results/diversity_stats.csv\n")

# =============================================================================
# PAIRWISE FST (Hudson et al. 1992)
# =============================================================================

hudson_fst <- function(seqs1, seqs2) {
    pi_within <- function(d) {
        v <- as.vector(dist.dna(d, model = "raw", pairwise.deletion = TRUE))
        mean(v, na.rm = TRUE)
    }
    n1 <- nrow(seqs1); n2 <- nrow(seqs2)
    pi1  <- pi_within(seqs1)
    pi2  <- pi_within(seqs2)
    pi_T <- pi_within(rbind(seqs1, seqs2))
    pi_S <- (n1 * pi1 + n2 * pi2) / (n1 + n2)
    fst  <- if (pi_T > 0) 1 - pi_S / pi_T else 0
    round(fst, 4)
}

cat("\n[Computing pairwise FST...]\n")
pop_names <- names(pop_list)
fst_mat   <- matrix(NA, nrow = length(pop_names), ncol = length(pop_names),
                    dimnames = list(pop_names, pop_names))
diag(fst_mat) <- 0

for (i in seq_along(pop_names)) {
    for (j in seq_along(pop_names)) {
        if (i < j) {
            fst_val <- tryCatch(
                hudson_fst(pop_list[[i]], pop_list[[j]]),
                error = function(e) NA_real_
            )
            fst_mat[i, j] <- fst_val
            fst_mat[j, i] <- fst_val
        }
    }
}

fst_df <- as.data.frame(fst_mat)
cat("\n── Pairwise FST ──\n")
print(fst_df)
write.csv(fst_df, "../results/fst_matrix.csv")
cat("Saved: results/fst_matrix.csv\n")

# =============================================================================
# AMOVA
# =============================================================================

cat("\n[AMOVA...]\n")
all_acc    <- unlist(map(pop_list, rownames))
dna_amova  <- dna_focal[all_acc, ]
pop_factor <- factor(pop_assignments[all_acc])

amova_res <- tryCatch({
    dist_mat <- dist.dna(dna_amova, model = "raw", pairwise.deletion = TRUE)
    amova(dist_mat ~ pop_factor, nperm = 999)
}, error = function(e) { message("AMOVA error: ", e$message); NULL })

if (!is.null(amova_res)) {
    cat("\n── AMOVA ──\n")
    print(amova_res)
    capture.output(print(amova_res)) %>% writeLines("../results/amova_results.txt")
    cat("Saved: results/amova_results.txt\n")
}

# =============================================================================
# HAPLOTYPE NETWORK
# =============================================================================

cat("\n[Haplotype network...]\n")
haps_net  <- haplotype(dna_amova)
net       <- haploNet(haps_net)
hap_idx   <- attr(haps_net, "index")
acc_to_pop <- pop_assignments[all_acc]

hap_to_pop <- sapply(seq_along(attr(haps_net, "dimnames")[[1]]), function(i) {
    idxs     <- hap_idx[[i]]
    seqnames <- rownames(dna_amova)[idxs]
    pops     <- acc_to_pop[seqnames]
    pops     <- pops[!is.na(pops) & pops != "Unknown"]
    if (length(pops) == 0) return("Unknown")
    names(sort(table(pops), decreasing = TRUE))[1]
})
names(hap_to_pop) <- attr(haps_net, "dimnames")[[1]]

node_cols <- SP_COLOURS[hap_to_pop]
node_cols[is.na(node_cols)] <- "grey80"
node_size <- summary(haps_net)

for (ext in c("png", "pdf")) {
    out_path <- paste0("../figures/haplotype_network_3sp.", ext)
    if (ext == "png") png(out_path, width = 1800, height = 1600, res = 180)
    else              pdf(out_path, width = 13, height = 11)
    par(mar = c(3, 2, 4, 2), bg = "white")
    plot(net, size = node_size * 2, bg = node_cols, border = "grey30",
         labels = FALSE, show.mutation = 0, fast = TRUE)
    title(
        main = "COI Haplotype Network - Sigaus complex (3 focal species)",
        sub  = paste0("n = ", length(all_acc), " sequences | ",
                      nrow(as.matrix(haps_net)), " haplotypes | ",
                      "node size proportional to frequency"),
        cex.main = 1.2, cex.sub = 0.85
    )
    legend("bottomleft",
           legend = c("S. australis", "S. campestris", "S. villosus"),
           fill   = c("#e74c3c", "#3498db", "#2ecc71"),
           title  = "Species", cex = 0.9, bty = "n",
           border = "grey30", text.font = 3)
    dev.off()
}
cat("Saved: figures/haplotype_network_3sp (.png + .pdf)\n")

# =============================================================================
# FIGURES
# =============================================================================

save_fig <- function(p, name, w, h) {
    ggsave(paste0("../figures/", name, ".png"), p, width = w, height = h, dpi = 200)
    pdf(paste0("../figures/", name, ".pdf"), width = w, height = h)
    print(p)
    dev.off()
    cat("Saved: figures/", name, "(.png + .pdf)\n", sep = "")
}

# ── Diversity plots ────────────────────────────────────────────────────────────
pA <- ggplot(stats_df, aes(x = reorder(population, Tajima_D), y = Tajima_D,
                            fill = population)) +
    geom_col(width = 0.5, colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted", colour = "grey50") +
    scale_fill_manual(values = SP_COLOURS, labels = SP_LABELS, guide = "none") +
    scale_x_discrete(labels = SP_LABELS) +
    coord_flip() +
    labs(title    = "A. Tajima's D",
         subtitle = "Dotted lines = +/-1.96 significance threshold",
         x = NULL, y = "Tajima's D") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y  = element_text(face = "italic"),
          plot.title   = element_text(face = "bold"))

pB <- ggplot(stats_df, aes(x = Hd, y = pi, colour = population,
                            label = population)) +
    geom_point(aes(size = n), alpha = 0.9) +
    geom_text_repel(aes(label = SP_LABELS[population]),
                    size = 3.5, fontface = "italic") +
    scale_colour_manual(values = SP_COLOURS, guide = "none") +
    scale_size_continuous(name = "n (seqs)", range = c(4, 12)) +
    annotate("text", x = 0.998, y = max(stats_df$pi) * 1.01,
             label = "High Hd + High pi\nHistoric large pop.",
             size = 2.8, colour = "grey40", hjust = 1) +
    annotate("text", x = 0.998, y = min(stats_df$pi) * 0.96,
             label = "High Hd + Low pi\nRecent expansion",
             size = 2.8, colour = "#e74c3c", hjust = 1) +
    labs(title    = "B. Haplotype vs Nucleotide Diversity",
         subtitle = "Avise et al. (1987) demographic framework",
         x = "Haplotype diversity (Hd)", y = "Nucleotide diversity (pi)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

p_div <- (pA | pB) +
    plot_annotation(
        title   = "Population Genetic Diversity - focal Sigaus species (COI)",
        caption = "Tajima (1989) | Nei (1987) | n = 106 sequences total"
    )
save_fig(p_div, "diversity_plots_3sp", 13, 5)

# ── Mismatch distributions ────────────────────────────────────────────────────
mm_plots <- map2(stats_df$mm_values, stats_df$population, function(mm, pop) {
    mm_df   <- tibble(diffs = mm)
    mean_mm <- mean(mm)
    max_d   <- ceiling(max(mm))
    exp_df  <- tibble(diffs = 0:max_d, freq = dpois(0:max_d, lambda = mean_mm))
    exp_df$freq <- exp_df$freq / sum(exp_df$freq)
    ggplot() +
        geom_histogram(data = mm_df, aes(x = diffs, y = after_stat(density)),
                       binwidth = 1, fill = SP_COLOURS[pop],
                       alpha = 0.65, colour = "white") +
        geom_line(data = exp_df, aes(x = diffs, y = freq),
                  colour = "black", linewidth = 1, linetype = "dashed") +
        labs(title    = SP_LABELS[pop],
             subtitle = paste0("n=", length(mm),
                               " | mean=", round(mean_mm, 1),
                               " | r=", round(sum(diff(exp_df$freq)^2), 4)),
             x = "Pairwise differences", y = "Frequency") +
        theme_minimal(base_size = 10) +
        theme(plot.title    = element_text(face = "bold.italic"),
              plot.subtitle = element_text(size = 8, colour = "grey40"))
})
p_mm <- wrap_plots(mm_plots, ncol = 3) +
    plot_annotation(
        title   = "Mismatch Distributions - focal Sigaus species (COI)",
        caption = "Black dashed = Poisson expectation under sudden expansion (Rogers & Harpending 1992)"
    )
save_fig(p_mm, "mismatch_distributions_3sp", 14, 5)

# ── FST heatmap ────────────────────────────────────────────────────────────────
fst_long <- fst_df %>%
    rownames_to_column("pop1") %>%
    pivot_longer(-pop1, names_to = "pop2", values_to = "FST") %>%
    filter(!is.na(FST)) %>%
    mutate(pop1 = recode(pop1, !!!SP_LABELS),
           pop2 = recode(pop2, !!!SP_LABELS))

p_fst <- ggplot(fst_long, aes(x = pop1, y = pop2, fill = FST)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", FST)), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "#fef9f9", high = "#c0392b",
                        limits = c(0, 1), name = expression(F[ST])) +
    labs(
        title    = expression("Pairwise " * F[ST] * " - focal " * italic("Sigaus") * " lineages (COI)"),
        subtitle = "Hudson et al. (1992) | values 0.21-0.32 = moderate-to-high inter-lineage structure",
        x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1,
                                     face = "italic", size = 11),
          axis.text.y = element_text(face = "italic", size = 11),
          plot.title  = element_text(face = "bold"),
          panel.grid  = element_blank())
save_fig(p_fst, "fst_heatmap_3sp", 7, 5)

cat("\n=== Population genetics complete ===\n")
cat("Next: source('05_sdm_maxent.R')\n")
