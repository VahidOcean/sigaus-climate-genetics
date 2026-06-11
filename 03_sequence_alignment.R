# =============================================================================
# 03_sequence_alignment.R
# =============================================================================
# Aligns raw COI sequences using MUSCLE and identifies unique haplotypes.
#
# Inputs:  data/sigaus_COI_raw.fasta
#          data/genbank_metadata.csv
# Outputs: data/sigaus_COI_aligned.fasta
#          data/haplotype_table.csv
#          data/haplotype_frequencies.csv
#          figures/alignment_summary.png / .pdf
#
# Runtime: ~2 minutes
# =============================================================================

library(tidyverse)
library(ape)
library(Biostrings)
library(muscle)
library(pegas)
library(patchwork)

cat("\n=== 03: Sequence Alignment + Haplotype Identification ===\n")

# ── Load raw sequences ────────────────────────────────────────────────────────
fasta_path <- "../data/sigaus_COI_raw.fasta"
if (!file.exists(fasta_path)) {
    stop(
        "FASTA not found: ", fasta_path, "\n",
        "Run 01_fetch_data.R first, or manually download sequences from:\n",
        "  https://www.ncbi.nlm.nih.gov/nuccore/\n",
        "  Query: Sigaus[ORGN] AND (COI[GENE] OR COX1[GENE])\n",
        "  Save as FASTA to data/sigaus_COI_raw.fasta"
    )
}

dna_raw <- readDNAStringSet(fasta_path)
names(dna_raw) <- sub("^(\\S+).*", "\\1", names(dna_raw))   # keep accession only
cat("Sequences loaded:", length(dna_raw), "\n")

# ── Length filter (standard COI barcode: 400–900 bp) ─────────────────────────
lens   <- width(dna_raw)
cat("Length range: min =", min(lens), "| max =", max(lens), "| median =", median(lens), "\n")
keep   <- lens >= 400 & lens <= 900
dna_filt <- dna_raw[keep]
cat("Retained (400–900 bp):", length(dna_filt), "/", length(dna_raw), "\n")

# ── MUSCLE alignment ──────────────────────────────────────────────────────────
cat("\n[Aligning with MUSCLE v3...]\n")
aln     <- muscle::muscle(dna_filt, quiet = FALSE)
aln_mat <- as.matrix(aln)
cat("Alignment: ", nrow(aln_mat), "seqs x", ncol(aln_mat), "columns\n")

# ── Trim alignment ────────────────────────────────────────────────────────────
# Remove columns with > 50% gaps
gap_freq  <- colMeans(aln_mat == "-")
aln_trim  <- aln_mat[, gap_freq < 0.5]
cat("After gap trimming:", ncol(aln_trim), "columns\n")

# Trim ragged ends (columns with < 5% coverage)
col_cov  <- colMeans(aln_trim != "-")
occupied <- which(col_cov >= 0.05)
aln_trim <- aln_trim[, min(occupied):max(occupied)]
cat("After end trimming:", ncol(aln_trim), "columns\n")

# Remove sequences with > 20% gaps after trimming
seq_gap  <- rowMeans(aln_trim == "-")
aln_trim <- aln_trim[seq_gap < 0.20, ]
cat("Sequences retained:", nrow(aln_trim), "\n")

# Write trimmed alignment
aln_strings <- DNAStringSet(apply(aln_trim, 1, paste, collapse = ""))
writeXStringSet(aln_strings, "../data/sigaus_COI_aligned.fasta")
cat("Saved: data/sigaus_COI_aligned.fasta\n")

# ── Haplotype identification ──────────────────────────────────────────────────
cat("\n[Identifying haplotypes...]\n")

dna_ape  <- read.dna("../data/sigaus_COI_aligned.fasta", format = "fasta")
haps     <- haplotype(dna_ape)
n_hap    <- nrow(as.matrix(haps))
n_seq    <- nrow(dna_ape)
cat("Sequences:", n_seq, "| Unique haplotypes:", n_hap, "\n")

# Accession -> haplotype mapping
ind2hap <- setNames(rep(NA_character_, n_seq), rownames(dna_ape))
for (h in rownames(as.matrix(haps))) {
    idxs         <- attr(haps, "index")[[which(rownames(as.matrix(haps)) == h)]]
    ind2hap[idxs] <- h
}

# ── Assign populations by organism name ──────────────────────────────────────
meta <- read_csv("../data/genbank_metadata.csv", show_col_types = FALSE)

assign_population <- function(organism) {
    case_when(
        str_detect(organism, regex("australis", ignore_case = TRUE)) ~ "S_australis_group",
        str_detect(organism, regex("campestris", ignore_case = TRUE)) ~ "S_campestris",
        str_detect(organism, regex("villosus", ignore_case = TRUE))   ~ "S_villosus",
        TRUE ~ "Unknown"
    )
}

hap_table <- tibble(accession = names(ind2hap), haplotype = unname(ind2hap)) %>%
    left_join(meta %>% select(accession, organism, seq_length), by = "accession") %>%
    mutate(population = assign_population(organism))

cat("\nHaplotypes per population:\n")
print(table(hap_table$population))

write_csv(hap_table, "../data/haplotype_table.csv")
cat("Saved: data/haplotype_table.csv\n")

# Haplotype frequencies
hap_freq <- hap_table %>%
    group_by(population, haplotype) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(population, desc(count))
write_csv(hap_freq, "../data/haplotype_frequencies.csv")
cat("Saved: data/haplotype_frequencies.csv\n")

# ── Figure: alignment summary ─────────────────────────────────────────────────
cat("\n[Figure] Alignment summary...\n")

len_df       <- tibble(length = sapply(as.list(dna_filt),
                                       function(x) length(as.vector(x))))
hap_count_df <- hap_table %>%
    group_by(haplotype) %>% summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
pop_hap_df   <- hap_table %>%
    group_by(population) %>%
    summarise(n_seq = n(), n_hap = n_distinct(haplotype), .groups = "drop")

pA <- ggplot(len_df, aes(x = length)) +
    geom_histogram(binwidth = 20, fill = "#2ecc71", colour = "white") +
    geom_vline(xintercept = c(400, 900), linetype = "dashed", colour = "red") +
    labs(title = "Raw sequence lengths",
         subtitle = "Red lines = 400 & 900 bp filter",
         x = "Length (bp)", y = "Count") +
    theme_minimal(base_size = 10)

pB <- ggplot(hap_count_df %>% slice_head(n = 20),
             aes(x = reorder(haplotype, n), y = n)) +
    geom_col(fill = "#3498db") + coord_flip() +
    labs(title = "Top 20 haplotype frequencies",
         x = "Haplotype", y = "Count") +
    theme_minimal(base_size = 10)

pC <- ggplot(pop_hap_df %>% filter(population != "Unknown"),
             aes(x = reorder(population, n_seq))) +
    geom_col(aes(y = n_seq), fill = "#e67e22", alpha = 0.7) +
    geom_point(aes(y = n_hap * max(n_seq) / max(n_hap)),
               colour = "#c0392b", size = 3) +
    scale_y_continuous(
        name = "No. sequences",
        sec.axis = sec_axis(
            ~ . * max(pop_hap_df$n_hap) / max(pop_hap_df$n_seq),
            name = "No. haplotypes"
        )
    ) +
    coord_flip() +
    labs(title = "Sequences (bars) & haplotypes (points) per group", x = NULL) +
    theme_minimal(base_size = 10)

p_aln <- (pA | pB) / pC +
    plot_annotation(
        title   = "Sigaus COI — Alignment Summary",
        caption = paste0("n = ", n_seq, " sequences | ",
                         n_hap, " unique haplotypes | ",
                         ncol(aln_trim), " alignment columns")
    )

ggsave("../figures/alignment_summary.png", p_aln,
       width = 12, height = 9, dpi = 200)
pdf("../figures/alignment_summary.pdf", width = 12, height = 9)
print(p_aln)
dev.off()
cat("Saved: figures/alignment_summary (.png + .pdf)\n")

cat("\n=== Alignment complete ===\n")
cat("Next: source('04_popgen_stats.R')\n")
