# =============================================================================
# 05_sdm_maxent.R
# =============================================================================
# Fits MaxEnt species distribution models for three focal Sigaus species
# and projects distributions under two CMIP6 climate change scenarios.
#
# Climate data:
#   Current:  WorldClim v2.1 baseline 1970-2000 (2.5 arc-min, ~5 km)
#   Future:   CMIP6 2081-2100 ensemble mean of three GCMs:
#             ACCESS-CM2, CMCC-ESM2, MRI-ESM2-0
#   Scenarios: SSP2-4.5 (+2 degrees C) and SSP5-8.5 (+4.4 degrees C)
#
# Bioclim variables used (selected for alpine ectotherms, low multicollinearity):
#   BIO1  = Annual Mean Temperature
#   BIO4  = Temperature Seasonality
#   BIO7  = Temperature Annual Range
#   BIO12 = Annual Precipitation
#   BIO15 = Precipitation Seasonality
#   BIO18 = Precipitation of Warmest Quarter
#
# MaxEnt settings:
#   Features: linear + hinge
#   Regularisation multiplier: 1.5
#   Replicates: 5-fold cross-validation
#   Output format: cloglog
#   Background: 10,000 random points within NZ extent
#
# Inputs:  data/gbif_clean.csv
# Outputs: results/sdm_evaluation.csv
#          results/variable_importance.csv
#          results/predicted_range_change.csv
#          results/predicted_range_change_summary.csv
#          figures/sdm_Sigaus_*.png / .pdf
#          figures/range_change_summary.png / .pdf
#
# Runtime: ~30 min first run (downloads ~2.6 GB climate data)
#          ~5 min subsequent runs (data cached)
# =============================================================================

library(tidyverse)
library(terra)
library(raster)
library(sf)
library(geodata)
library(dismo)
library(rnaturalearth)
library(patchwork)
library(viridis)
library(rJava)

cat("\n=== 05: Species Distribution Modelling (MaxEnt) ===\n")

# ── Configuration ─────────────────────────────────────────────────────────────
NZ_EXT       <- ext(166, 174.5, -47.5, -34)
SI_LAT_RANGE <- c(-47.5, -40)
BIOCLIM_VARS <- c(1, 4, 7, 12, 15, 18)
RES          <- "2.5"
FUTURE_YEAR  <- "2081-2100"
GCM_MODELS   <- c("ACCESS-CM2", "CMCC-ESM2", "MRI-ESM2-0")
SCENARIOS    <- c("ssp245", "ssp585")
THRESHOLD    <- 0.5          # cloglog threshold for binary maps

dir.create("../data/climate", showWarnings = FALSE, recursive = TRUE)

# Check MaxEnt jar
maxent_jar <- file.path(system.file(package = "dismo"), "java", "maxent.jar")
if (!file.exists(maxent_jar)) {
    stop(
        "maxent.jar not found at: ", maxent_jar, "\n",
        "Run shell/00_setup.sh to install it, or download from:\n",
        "  https://biodiversityinformatics.amnh.org/open_source/maxent/\n",
        "Place in: ", dirname(maxent_jar)
    )
}
.jinit()

# =============================================================================
# PART 1: Download climate data
# =============================================================================
cat("\n[1] Downloading WorldClim climate data...\n")

wc_all     <- worldclim_global(var = "bio", res = RES, path = "../data/climate")
wc_current <- wc_all[[BIOCLIM_VARS]]
wc_current <- crop(wc_current, NZ_EXT)
cat("  Current climate: ", nlyr(wc_current), "layers loaded\n")

cat("\n[2] Downloading future climate layers (CMIP6)...\n")
future_stacks <- list()

for (scenario in SCENARIOS) {
    gcm_layers <- list()
    for (gcm in GCM_MODELS) {
        cat("  Downloading:", gcm, scenario, FUTURE_YEAR, "\n")
        fut <- tryCatch({
            cmip6_world(
                model = gcm, ssp = sub("ssp", "", scenario),
                time  = FUTURE_YEAR, var = "bioc",
                res   = RES, path = "../data/climate"
            )[[BIOCLIM_VARS]]
        }, error = function(e) {
            message("  WARNING: Could not download ", gcm, " ", scenario)
            NULL
        })
        if (!is.null(fut)) gcm_layers[[gcm]] <- crop(fut, NZ_EXT)
        Sys.sleep(1)
    }
    if (length(gcm_layers) > 0) {
        ensemble <- rast(lapply(seq_len(nlyr(gcm_layers[[1]])), function(i) {
            mean(rast(lapply(gcm_layers, `[[`, i)))
        }))
        names(ensemble) <- names(wc_current)
        future_stacks[[scenario]] <- ensemble
        cat("  Ensemble ready for", scenario, "(", length(gcm_layers), "GCMs)\n")
    }
}

# Convert to RasterStack for dismo compatibility
wc_raster     <- raster::stack(wc_current)
future_rasters <- map(future_stacks, raster::stack)

# =============================================================================
# PART 2: Prepare occurrence data
# =============================================================================
cat("\n[3] Preparing occurrence data...\n")

occ <- read_csv("../data/gbif_clean.csv", show_col_types = FALSE)
names(occ)[names(occ) == "decimallatitude"]  <- "lat"
names(occ)[names(occ) == "decimallongitude"] <- "lon"
occ <- occ %>% filter(!is.na(lat), !is.na(lon))

species_counts <- dplyr::count(occ, species) %>% arrange(desc(n))
cat("Records per species:\n"); print(species_counts)

model_taxa <- species_counts %>% filter(n >= 10) %>% pull(species)
cat("Taxa to model:", paste(model_taxa, collapse = ", "), "\n")

# Background points
set.seed(42)
bg_pts <- terra::spatSample(wc_current[[1]], size = 10000,
                             method = "random", na.rm = TRUE, xy = TRUE)[, c("x","y")]
names(bg_pts) <- c("lon", "lat")
cat("Background points:", nrow(bg_pts), "\n")

# =============================================================================
# PART 3: Run MaxEnt models
# =============================================================================
cat("\n[4] Running MaxEnt models...\n")

sdm_results  <- list()
eval_metrics <- list()
var_imp_list <- list()

for (sp in model_taxa) {
    cat("\n  Species:", sp, "\n")

    sp_occ  <- occ %>% filter(species == sp) %>%
        dplyr::select(lon, lat) %>% as.data.frame()
    occ_env <- as.data.frame(raster::extract(wc_raster, sp_occ))
    bg_env  <- as.data.frame(raster::extract(wc_raster, bg_pts))

    occ_clean <- sp_occ[complete.cases(occ_env), ]
    bg_clean  <- bg_pts[complete.cases(bg_env), ]
    cat("  Clean occurrences:", nrow(occ_clean), "\n")
    if (nrow(occ_clean) < 5) { cat("  Too few — skipping\n"); next }

    out_dir <- paste0("../results/maxent_", gsub(" ", "_", sp))
    dir.create(out_dir, showWarnings = FALSE)

    mx <- tryCatch(
        dismo::maxent(
            x    = wc_raster,
            p    = occ_clean,
            a    = bg_clean,
            path = out_dir,
            args = c(
                "linear=true", "quadratic=false", "product=false",
                "threshold=false", "hinge=true",
                "betamultiplier=1.5",
                "replicates=5", "replicatetype=crossvalidate",
                "outputformat=cloglog",
                "responsecurves=true", "jackknife=true"
            )
        ),
        error = function(e) { message("  MaxEnt error: ", e$message); NULL }
    )
    if (is.null(mx)) next

    # Model evaluation
    eval_file <- list.files(out_dir, pattern = "maxentResults.csv", full.names = TRUE)
    if (length(eval_file) > 0) {
        ec       <- read.csv(eval_file[1])
        tr_auc   <- mean(ec$Training.AUC, na.rm = TRUE)
        te_auc   <- mean(ec$Test.AUC,     na.rm = TRUE)
        cat("  Train AUC:", round(tr_auc, 3), "| Test AUC:", round(te_auc, 3), "\n")
        eval_metrics[[sp]] <- tibble(
            species      = sp,
            n_occ        = nrow(occ_clean),
            training_AUC = round(tr_auc, 3),
            test_AUC     = round(te_auc, 3)
        )
    }

    # Variable importance
    var_imp_list[[sp]] <- tryCatch({
        vi <- var.importance(mx); vi$species <- sp; vi
    }, error = function(e) NULL)

    # Predictions
    cat("  Predicting current and future distributions...\n")
    sdm_results[[sp]] <- list(
        model   = mx,
        current = predict(mx, wc_raster),
        future  = map(future_rasters, function(fut) {
            tryCatch(predict(mx, fut), error = function(e) NULL)
        })
    )
    cat("  Done:", sp, "\n")
}

# Save evaluation and variable importance
if (length(eval_metrics) > 0) {
    bind_rows(eval_metrics) %>% write_csv("../results/sdm_evaluation.csv")
    cat("\nSaved: results/sdm_evaluation.csv\n")
}
if (length(var_imp_list) > 0) {
    bind_rows(var_imp_list) %>% write_csv("../results/variable_importance.csv")
    cat("Saved: results/variable_importance.csv\n")
}

# =============================================================================
# PART 4: Range change calculations
# =============================================================================
cat("\n[5] Calculating range area changes...\n")

range_change <- map_dfr(names(sdm_results), function(sp) {
    res     <- sdm_results[[sp]]
    cur_bin <- res$current >= THRESHOLD
    area_cur <- cellStats(cur_bin, sum)
    map_dfr(names(res$future), function(scen) {
        fut <- res$future[[scen]]
        if (is.null(fut)) return(NULL)
        fut_bin <- fut >= THRESHOLD
        tibble(
            species      = sp,
            scenario     = scen,
            area_current = area_cur,
            area_future  = cellStats(fut_bin, sum),
            area_gained  = cellStats((fut_bin == 1) & (cur_bin == 0), sum),
            area_lost    = cellStats((cur_bin == 1) & (fut_bin == 0), sum),
            area_stable  = cellStats((cur_bin == 1) & (fut_bin == 1), sum),
            pct_change   = round((cellStats(fut_bin, sum) - area_cur) / area_cur * 100, 1)
        )
    })
})

write_csv(range_change, "../results/predicted_range_change.csv")

# Average across cross-validation replicates
range_summary <- range_change %>%
    group_by(species, scenario) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
write_csv(range_summary, "../results/predicted_range_change_summary.csv")

cat("\n── Mean range change ──\n")
print(range_summary %>% dplyr::select(species, scenario, pct_change))
cat("Saved: results/predicted_range_change_summary.csv\n")

# =============================================================================
# PART 5: Figures
# =============================================================================
cat("\n[6] Generating SDM figures...\n")

nz_coast <- ne_countries(country = "New Zealand", scale = "medium", returnclass = "sf")

save_fig <- function(p, name, w, h) {
    ggsave(paste0("../figures/", name, ".png"), p, width = w, height = h, dpi = 200)
    pdf(paste0("../figures/", name, ".pdf"), width = w, height = h)
    print(p); dev.off()
    cat("Saved: figures/", name, "(.png + .pdf)\n", sep = "")
}

# SDM maps per species
for (sp in names(sdm_results)) {
    res      <- sdm_results[[sp]]
    cur_mean <- mean(res$current)
    cur_df   <- raster::as.data.frame(cur_mean, xy = TRUE)
    names(cur_df) <- c("lon", "lat", "suitability")
    cur_df <- cur_df[!is.na(cur_df$suitability) &
                     cur_df$lat <= SI_LAT_RANGE[2] &
                     cur_df$lat >= SI_LAT_RANGE[1], ]
    occ_sp <- occ[occ$species == sp & occ$lat <= SI_LAT_RANGE[2], ]

    p_cur <- ggplot() +
        geom_raster(data = cur_df, aes(x = lon, y = lat, fill = suitability)) +
        geom_sf(data = nz_coast, fill = NA, colour = "black", linewidth = 0.5) +
        geom_point(data = occ_sp, aes(x = lon, y = lat),
                   shape = 3, size = 1.2, colour = "white", alpha = 0.7) +
        scale_fill_viridis_c(option = "inferno", name = "Suitability", limits = c(0,1)) +
        coord_sf(xlim = c(166, 174.5), ylim = SI_LAT_RANGE) +
        labs(title    = "Current habitat suitability",
             subtitle = paste0(sp, " | baseline 1970-2000")) +
        theme_minimal(base_size = 10) +
        theme(axis.title    = element_blank(),
              plot.title    = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))

    fut_plots <- imap(res$future, function(fut_rast, scen) {
        if (is.null(fut_rast)) return(NULL)
        diff_df   <- raster::as.data.frame(mean(fut_rast) - cur_mean, xy = TRUE)
        names(diff_df) <- c("lon", "lat", "change")
        diff_df <- diff_df[!is.na(diff_df$change) &
                           diff_df$lat <= SI_LAT_RANGE[2] &
                           diff_df$lat >= SI_LAT_RANGE[1], ]
        scen_lbl <- ifelse(scen == "ssp245", "SSP2-4.5 (+2 C)", "SSP5-8.5 (+4.4 C)")
        pct <- range_summary %>%
            filter(species == sp, scenario == scen) %>%
            pull(pct_change) %>% mean() %>% round(1)
        ggplot() +
            geom_raster(data = diff_df, aes(x = lon, y = lat, fill = change)) +
            geom_sf(data = nz_coast, fill = NA, colour = "black", linewidth = 0.5) +
            scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#27ae60",
                                 midpoint = 0, limits = c(-1, 1),
                                 name = "Delta Suit.") +
            coord_sf(xlim = c(166, 174.5), ylim = SI_LAT_RANGE) +
            labs(title    = paste("Change under", scen_lbl),
                 subtitle = paste0("Mean range change: ", pct, "%")) +
            theme_minimal(base_size = 10) +
            theme(axis.title    = element_blank(),
                  plot.title    = element_text(face = "bold"),
                  plot.subtitle = element_text(colour = "#c0392b", face = "bold"))
    }) %>% keep(~!is.null(.x))

    p_sdm <- (p_cur | fut_plots[[1]] | fut_plots[[2]]) +
        plot_annotation(
            title   = paste("SDM -", sp),
            caption = paste0("MaxEnt | WorldClim 2.1 | CMIP6 ensemble (3 GCMs) | ",
                             "5-fold CV | Test AUC = ",
                             eval_metrics[[sp]]$test_AUC),
            theme   = theme(plot.title = element_text(face = "bold.italic", size = 13))
        )
    save_fig(p_sdm, paste0("sdm_", gsub(" ", "_", sp)), 15, 6)
}

# Range change summary bar chart
p_range <- ggplot(range_summary,
                  aes(x = reorder(species, pct_change),
                      y = pct_change, fill = scenario)) +
    geom_col(position = position_dodge(0.7), width = 0.6,
             colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(
        values = c("ssp245" = "#f39c12", "ssp585" = "#c0392b"),
        labels = c("ssp245" = "SSP2-4.5 (+2 C)", "ssp585" = "SSP5-8.5 (+4.4 C)"),
        name   = "Climate scenario"
    ) +
    coord_flip() +
    labs(
        title    = "Predicted range area change - NZ alpine grasshoppers",
        subtitle = "2081-2100 vs baseline 1970-2000 | MaxEnt | CMIP6 ensemble (3 GCMs)",
        x        = NULL,
        y        = "Mean % change in suitable area",
        caption  = "Threshold = 0.5 cloglog | 5-fold cross-validation mean"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title  = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic"),
          legend.position = "right")
save_fig(p_range, "range_change_summary", 10, 5)

cat("\n=== SDM complete ===\n")
cat("Next: source('06_integrate_results.R')\n")
