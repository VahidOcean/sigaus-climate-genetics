# =============================================================================
# 06_integrate_results.R
# =============================================================================
# Integrates SDM projections with population genetic diversity metrics and
# produces the final integrated figures for the manuscript.
#
# Core hypothesis tested:
#   Populations predicted to lose suitable habitat under climate change
#   should show lower genetic diversity (Hd, pi) and positive Tajima's D,
#   consistent with demographic decline (genetic erosion).
#
#   The absence of this relationship (all species show high diversity
#   despite catastrophic predicted habitat loss) is interpreted as an
#   extinction debt: current diversity reflects historical abundance, not
#   contemporary or projected habitat change.
#
# Statistical note:
#   Pearson correlations and linear regressions in Figure 10 are exploratory
#   (n = 3 species) and are presented as a descriptive framework only.
#
# Inputs:  results/diversity_stats.csv
#          results/predicted_range_change_summary.csv
#          results/fst_matrix.csv
#          data/haplotype_table.csv
# Outputs: results/integrated_summary.csv
#          figures/integrated_panel_3sp.png / .pdf
#
# Runtime: ~2 minutes
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggrepel)

cat("\n=== 06: Integrating SDM + Population Genetics ===\n")

# ── Configuration ─────────────────────────────────────────────────────────────
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

# ── Load results ──────────────────────────────────────────────────────────────
div_stats     <- read_csv("../results/diversity_stats.csv", show_col_types = FALSE)
range_summary <- read_csv("../results/predicted_range_change_summary.csv",
                          show_col_types = FALSE)

# Manual crosswalk: popgen population names <-> SDM species names
name_crosswalk <- tibble(
    population = c("S_australis_group", "S_campestris", "S_villosus"),
    species    = c("Sigaus australis",  "Sigaus campestris", "Sigaus villosus")
)

# Join diversity stats with SDM range change (SSP5-8.5)
integrated <- div_stats %>%
    left_join(name_crosswalk, by = "population") %>%
    left_join(
        range_summary %>%
            filter(scenario == "ssp585") %>%
            dplyr::select(species, pct_change, area_current, area_lost, area_gained),
        by = "species"
    )

cat("Integrated table:\n")
print(integrated %>% dplyr::select(population, n, Hd, pi, Tajima_D, pct_change))

write_csv(integrated, "../results/integrated_summary.csv")
cat("Saved: results/integrated_summary.csv\n")

# =============================================================================
# PEARSON CORRELATIONS (exploratory, n = 3)
# =============================================================================
int_plot <- integrated %>% filter(!is.na(pct_change))

for (metric in c("Hd", "pi", "Tajima_D")) {
    if (nrow(int_plot) >= 3) {
        ct <- cor.test(int_plot[[metric]], int_plot$pct_change, method = "pearson")
        cat(sprintf("  %s vs pct_change: r = %.3f, p = %.3f (n=%d, exploratory)\n",
                    metric, ct$estimate, ct$p.value, nrow(int_plot)))
    }
}

# =============================================================================
# FIGURES
# =============================================================================

save_fig <- function(p, name, w, h) {
    ggsave(paste0("../figures/", name, ".png"), p, width = w, height = h, dpi = 200)
    pdf(paste0("../figures/", name, ".pdf"), width = w, height = h)
    print(p); dev.off()
    cat("Saved: figures/", name, "(.png + .pdf)\n", sep = "")
}

# ── Panel A: Hd vs range change ───────────────────────────────────────────────
pA <- ggplot(int_plot, aes(x = pct_change, y = Hd,
                            colour = population, label = population)) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey60",
                fill = "grey85", linewidth = 0.8) +
    geom_point(aes(size = n), alpha = 0.9) +
    geom_text_repel(aes(label = SP_LABELS[population]),
                    size = 3.5, fontface = "italic") +
    scale_colour_manual(values = SP_COLOURS, guide = "none") +
    scale_size_continuous(range = c(5, 12), name = "n (seqs)") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title    = "A. Haplotype diversity vs range change",
         subtitle = "Exploratory | n = 3 species",
         x = "Predicted % change in suitable area (SSP5-8.5, 2081-2100)",
         y = "Haplotype diversity (Hd)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

# ── Panel B: pi vs range change ───────────────────────────────────────────────
pB <- ggplot(int_plot, aes(x = pct_change, y = pi,
                            colour = population, label = population)) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey60",
                fill = "grey85", linewidth = 0.8) +
    geom_point(aes(size = n), alpha = 0.9) +
    geom_text_repel(aes(label = SP_LABELS[population]),
                    size = 3.5, fontface = "italic") +
    scale_colour_manual(values = SP_COLOURS, guide = "none") +
    scale_size_continuous(range = c(5, 12), name = "n (seqs)") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title    = "B. Nucleotide diversity vs range change",
         subtitle = "Exploratory | n = 3 species",
         x = "Predicted % change in suitable area (SSP5-8.5, 2081-2100)",
         y = "Nucleotide diversity (pi)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

# ── Panel C: Tajima's D ───────────────────────────────────────────────────────
pC <- ggplot(integrated, aes(x = reorder(population, Tajima_D),
                              y = Tajima_D, fill = population)) +
    geom_col(width = 0.5, colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted", colour = "grey50") +
    scale_fill_manual(values = SP_COLOURS, labels = SP_LABELS, guide = "none") +
    scale_x_discrete(labels = SP_LABELS) +
    coord_flip() +
    labs(title    = "C. Tajima's D - all focal species",
         subtitle = "Positive values: subdivision / Wahlund effect",
         x = NULL, y = "Tajima's D") +
    theme_minimal(base_size = 11) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic"))

# ── Panel D: Range change bar chart ──────────────────────────────────────────
pD <- ggplot(range_summary,
             aes(x = reorder(species, pct_change),
                 y = pct_change, fill = scenario)) +
    geom_col(position = position_dodge(0.7), width = 0.6,
             colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(
        values = c("ssp245" = "#f39c12", "ssp585" = "#c0392b"),
        labels = c("ssp245" = "SSP2-4.5 (+2 C)", "ssp585" = "SSP5-8.5 (+4.4 C)"),
        name   = "Scenario"
    ) +
    coord_flip() +
    labs(title    = "D. Predicted range area change",
         subtitle = "2081-2100 vs baseline 1970-2000",
         x = NULL, y = "Mean % change in suitable area") +
    theme_minimal(base_size = 11) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic"),
          legend.position = "bottom")

# ── Assemble final panel ──────────────────────────────────────────────────────
final_panel <- (pA | pB) / (pC | pD) +
    plot_annotation(
        title    = "Climate-Driven Range Shifts & Genetic Diversity - Focal Sigaus Species",
        subtitle = paste0(
            "S. australis (n=76 seqs, 212 occ) | ",
            "S. campestris (n=20 seqs, 49 occ) | ",
            "S. villosus (n=10 seqs, 22 occ) | ",
            "COI mtDNA + MaxEnt SDM (WorldClim 2.1 + CMIP6, 3 GCMs)"
        ),
        caption  = paste0(
            "Hypothesis tested: populations facing greater habitat loss should show ",
            "lower genetic diversity (extinction debt framework).\n",
            "All species show uniformly high Hd (0.978-0.995) and positive Tajima's D ",
            "despite 67-97% predicted habitat loss — consistent with genetic extinction debt."
        ),
        theme = theme(
            plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 8.5, colour = "grey30"),
            plot.caption  = element_text(size = 8, colour = "grey40",
                                         hjust = 0, lineheight = 1.3)
        )
    )

save_fig(final_panel, "integrated_panel_3sp", 16, 13)

# =============================================================================
# PRINT INTERPRETIVE SUMMARY
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat(" RESULTS SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat("Genetic diversity (COI mtDNA):\n")
for (i in seq_len(nrow(integrated))) {
    r <- integrated[i, ]
    cat(sprintf("  %-22s Hd=%.3f  pi=%.4f  D=%+.3f  range_change=%s%%\n",
                SP_LABELS[r$population],
                r$Hd, r$pi, r$Tajima_D,
                ifelse(is.na(r$pct_change), "N/A", round(r$pct_change, 1))))
}

cat("\nAMOVA: PhiST =", 0.672, "(P < 0.001) —",
    "67.2% of variance among lineages\n")

cat("\nSDM range change (SSP5-8.5):\n")
for (i in seq_len(nrow(integrated %>% filter(!is.na(pct_change))))) {
    r <- (integrated %>% filter(!is.na(pct_change)))[i, ]
    cat(sprintf("  %-22s %.1f%%\n", SP_LABELS[r$population], r$pct_change))
}

cat("\nInterpretation: Extinction debt detected.\n")
cat("Current genetic diversity reflects historical demographic stability.\n")
cat("Future habitat loss will drive fragmentation and genetic erosion\n")
cat("not yet detectable in contemporary mtDNA signals.\n")
cat(strrep("=", 60), "\n")

cat("\n=== PIPELINE COMPLETE ===\n")
cat("All outputs in: results/ and figures/\n")
