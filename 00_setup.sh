#!/usr/bin/env bash
# =============================================================================
# 00_setup.sh  —  One-time environment setup
# =============================================================================
# Run this script ONCE from Terminal before opening RStudio.
# Installs all system libraries and R packages required by the pipeline.
#
# Usage:
#   chmod +x shell/00_setup.sh
#   ./shell/00_setup.sh
#
# Tested on: macOS 14+ (Apple Silicon and Intel), Ubuntu 22.04+
# Requires:  Homebrew (macOS), R >= 4.5, Java >= 11
# =============================================================================

set -e
echo "============================================================"
echo " Sigaus Range-Shift Pipeline — Environment Setup"
echo "============================================================"

# ── 1. Homebrew (macOS only) ──────────────────────────────────────────────
if [[ "$OSTYPE" == "darwin"* ]]; then
    if ! command -v brew &>/dev/null; then
        echo "[1/4] Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    else
        echo "[1/4] Homebrew already installed — skipping."
    fi

    echo "[2/4] Installing system libraries (gdal, geos, proj, udunits, mafft)..."
    brew install gdal proj geos udunits mafft openjdk@17

    # Link Java for R CMD javareconf
    sudo ln -sfn /opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk \
        /Library/Java/JavaVirtualMachines/openjdk-17.jdk 2>/dev/null || true
    sudo R CMD javareconf 2>/dev/null || true

else
    echo "[1/4] Linux detected — using apt..."
    sudo apt-get update -qq
    sudo apt-get install -y \
        libgdal-dev libproj-dev libgeos-dev libudunits2-dev \
        mafft default-jdk r-base
    echo "[2/4] System libraries installed."
fi

# ── 3. R packages ─────────────────────────────────────────────────────────
echo "[3/4] Installing R packages (this may take 10–15 min first time)..."

Rscript - <<'REOF'

install_if_missing <- function(pkgs, repos = "https://cloud.r-project.org") {
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]
    if (length(missing) > 0) {
        message("Installing: ", paste(missing, collapse = ", "))
        install.packages(missing, repos = repos, dependencies = TRUE)
    } else {
        message("All CRAN packages already installed.")
    }
}

# CRAN packages
cran_pkgs <- c(
    # Data wrangling
    "tidyverse", "readr", "dplyr", "stringr", "lubridate",
    # Spatial / SDM
    "terra", "sf", "raster", "geodata", "dismo",
    "ENMeval", "CoordinateCleaner", "rgbif",
    "rnaturalearth", "rnaturalearthdata",
    # Population genetics
    "ape", "pegas", "adegenet", "hierfstat", "mmod", "seqinr",
    # Bioinformatics
    "rentrez",
    # Visualisation
    "ggplot2", "patchwork", "ggrepel", "viridis", "scales"
)
install_if_missing(cran_pkgs)

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("Biostrings", "muscle"), ask = FALSE, update = FALSE)

# rJava (needed for MaxEnt via dismo)
if (!"rJava" %in% rownames(installed.packages())) {
    install.packages("rJava", repos = "https://cloud.r-project.org")
}

message("\n All R packages installed successfully.")
REOF

# ── 4. MaxEnt jar ─────────────────────────────────────────────────────────
echo "[4/4] Checking MaxEnt jar..."

MAXENT_DIR="$(Rscript -e 'cat(system.file(package="dismo"))')/java"
mkdir -p "$MAXENT_DIR"

if [ ! -f "$MAXENT_DIR/maxent.jar" ]; then
    echo "  Downloading maxent.jar..."
    curl -L \
        "https://github.com/mrmaxent/Maxent/blob/master/ArchivedReleases/3.4.1/maxent.jar?raw=true" \
        -o "$MAXENT_DIR/maxent.jar"
    echo "  maxent.jar installed at: $MAXENT_DIR/maxent.jar"
else
    echo "  maxent.jar already present — skipping."
fi

# ── Done ──────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " Setup complete!"
echo ""
echo " Open RStudio and run scripts in order:"
echo "   setwd('path/to/sigaus_range_genetics/R')"
echo "   source('01_fetch_data.R')"
echo "   source('02_clean_occurrences.R')"
echo "   source('03_sequence_alignment.R')"
echo "   source('04_popgen_stats.R')"
echo "   source('05_sdm_maxent.R')"
echo "   source('06_integrate_results.R')"
echo "============================================================"
