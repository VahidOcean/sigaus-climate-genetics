# =============================================================================
# 02_clean_occurrences.R
# =============================================================================
# Applies a standardised coordinate cleaning protocol to raw GBIF records.
#
# Cleaning steps:
#   1. CoordinateCleaner automated flags (ocean, centroid, institution, zeros)
#   2. Coordinate uncertainty filter (> 10 km removed)
#   3. Temporal filter (pre-1950 records removed)
#   4. Spatial thinning (1 record per 5-km grid cell)
#
# Inputs:  data/gbif_raw.csv
# Outputs: data/gbif_clean.csv
#          figures/occurrences_map_3sp.png
#          figures/occurrences_map_3sp.pdf
#
# Runtime: ~1 minute
# =============================================================================

library(tidyverse)
library(CoordinateCleaner)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

cat("\n=== 02: Cleaning GBIF occurrences ===\n")

# ── Load data ─────────────────────────────────────────────────────────────────
gbif_raw <- read_csv("../data/gbif_raw.csv", show_col_types = FALSE)
cat("Raw records:", nrow(gbif_raw), "\n")

# Standardise column names for CoordinateCleaner
occ <- gbif_raw %>%
    rename(
        decimallongitude = decimalLongitude,
        decimallatitude  = decimalLatitude
    ) %>%
    filter(!is.na(decimallongitude), !is.na(decimallatitude))

cat("Records with coordinates:", nrow(occ), "\n")

# ── Step 1: CoordinateCleaner ─────────────────────────────────────────────────
cat("\n[Step 1] CoordinateCleaner automated flags...\n")

cc_flags <- clean_coordinates(
    x       = occ,
    lon     = "decimallongitude",
    lat     = "decimallatitude",
    species = "species",
    tests   = c("capitals", "centroids", "equal", "gbif",
                "institutions", "seas", "zeros"),
    value   = "spatialvalid"
)
cat("  Flagged:", sum(!cc_flags$.summary), "records\n")
cat("  Passing:", sum(cc_flags$.summary), "records\n")
occ_cc <- occ[cc_flags$.summary, ]

# ── Step 2: Coordinate uncertainty ───────────────────────────────────────────
cat("\n[Step 2] Coordinate uncertainty filter (> 10 km)...\n")
n_before <- nrow(occ_cc)
if ("coordinateUncertaintyInMeters" %in% names(occ_cc)) {
    occ_cc <- occ_cc %>%
        filter(is.na(coordinateUncertaintyInMeters) |
               coordinateUncertaintyInMeters <= 10000)
}
cat("  Removed:", n_before - nrow(occ_cc), "records\n")

# ── Step 3: Temporal filter ───────────────────────────────────────────────────
cat("\n[Step 3] Temporal filter (>= 1950)...\n")
n_before <- nrow(occ_cc)
occ_cc   <- occ_cc %>% filter(is.na(year) | year >= 1950)
cat("  Removed:", n_before - nrow(occ_cc), "records\n")

# ── Step 4: Spatial thinning ──────────────────────────────────────────────────
cat("\n[Step 4] Spatial thinning (1 record per ~5-km cell)...\n")
occ_thin <- occ_cc %>%
    mutate(
        cell_lat = round(decimallatitude,  2),
        cell_lon = round(decimallongitude, 2)
    ) %>%
    group_by(species, cell_lat, cell_lon) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-cell_lat, -cell_lon) %>%
    rename(lat = decimallatitude, lon = decimallongitude)

cat("  Before thinning:", nrow(occ_cc), "\n")
cat("  After thinning: ", nrow(occ_thin), "\n")

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\nClean records per species:\n")
print(table(occ_thin$species))

write_csv(occ_thin, "../data/gbif_clean.csv")
cat("\nSaved: data/gbif_clean.csv\n")

# ── Figure: occurrence map ────────────────────────────────────────────────────
cat("\n[Figure] Generating occurrence map...\n")

nz <- ne_countries(country = "New Zealand", scale = "medium", returnclass = "sf")

occ_sf <- st_as_sf(occ_thin, coords = c("lon", "lat"), crs = 4326)

pal_3 <- c(
    "Sigaus australis"  = "#e74c3c",
    "Sigaus campestris" = "#3498db",
    "Sigaus villosus"   = "#2ecc71"
)

p_occ <- ggplot() +
    geom_sf(data = nz, fill = "grey92", colour = "grey60", linewidth = 0.3) +
    geom_sf(data = occ_sf, aes(colour = species),
            size = 1.8, alpha = 0.8, shape = 16) +
    scale_colour_manual(values = pal_3, name = "Species") +
    coord_sf(xlim = c(166, 174.5), ylim = c(-47.5, -40)) +
    labs(
        title    = expression("Focal " * italic("Sigaus") * " species — Clean GBIF Occurrences"),
        subtitle = paste0("n = ", nrow(occ_thin),
                          " records after coordinate cleaning & 5-km spatial thinning"),
        caption  = "Source: GBIF | CoordinateCleaner (Zizka et al. 2019)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.text      = element_text(face = "italic", size = 9),
        plot.title       = element_text(face = "bold"),
        panel.grid.major = element_line(colour = "grey85", linewidth = 0.3)
    )

ggsave("../figures/occurrences_map_3sp.png", p_occ,
       width = 8, height = 9, dpi = 200)
pdf("../figures/occurrences_map_3sp.pdf", width = 8, height = 9)
print(p_occ)
dev.off()
cat("Saved: figures/occurrences_map_3sp (.png + .pdf)\n")

cat("\n=== Cleaning complete ===\n")
cat("Next: source('03_sequence_alignment.R')\n")
