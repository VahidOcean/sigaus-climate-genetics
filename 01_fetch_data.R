# =============================================================================
# 01_fetch_data.R
# =============================================================================
# Downloads occurrence records (GBIF) and COI sequences (NCBI GenBank) for
# three focal Sigaus grasshopper species.
#
# Inputs:  none (downloads from public APIs)
# Outputs: data/gbif_raw.csv
#          data/genbank_metadata.csv
#          data/sigaus_COI_raw.fasta
#
# Runtime: ~5 minutes (depends on internet speed)
# =============================================================================

library(tidyverse)
library(rgbif)
library(rentrez)

# ── Configuration ─────────────────────────────────────────────────────────────
Entrez.email <- "your_email@example.com"   # <-- set your email for NCBI

FOCAL_SPECIES <- c(
    "Sigaus australis",
    "Sigaus campestris",
    "Sigaus villosus"
)

# NCBI search queries per marker
GENBANK_QUERIES <- list(
    COI = '(Sigaus[ORGN]) AND (COI[GENE] OR COX1[GENE] OR "cytochrome oxidase"[TITL])'
)

# ── Create output directories ──────────────────────────────────────────────────
dir.create("../data",    showWarnings = FALSE, recursive = TRUE)
dir.create("../results", showWarnings = FALSE)
dir.create("../figures", showWarnings = FALSE)

# =============================================================================
# PART 1: GBIF occurrences
# =============================================================================
cat("\n=== PART 1: Downloading GBIF occurrences ===\n")

fetch_gbif_taxon <- function(species_name, limit = 500) {
    cat("  Fetching:", species_name, "\n")
    name_result <- name_backbone(name = species_name, rank = "species")
    if (is.null(name_result$usageKey)) {
        warning("No GBIF backbone key for: ", species_name)
        return(NULL)
    }
    occ <- occ_search(
        taxonKey           = name_result$usageKey,
        hasCoordinate      = TRUE,
        hasGeospatialIssue = FALSE,
        country            = "NZ",
        limit              = limit,
        fields             = c(
            "key", "species", "decimalLatitude", "decimalLongitude",
            "coordinateUncertaintyInMeters", "elevation",
            "year", "month", "basisOfRecord",
            "institutionCode", "collectionCode", "catalogNumber",
            "locality", "datasetName", "stateProvince"
        )
    )
    if (is.null(occ$data) || nrow(occ$data) == 0) return(NULL)
    cat("    Retrieved:", nrow(occ$data), "records\n")
    occ$data
}

gbif_list <- map(FOCAL_SPECIES, function(sp) {
    result <- tryCatch(fetch_gbif_taxon(sp), error = function(e) {
        message("  ERROR: ", sp, " — ", e$message)
        NULL
    })
    Sys.sleep(0.5)
    result
})

gbif_raw <- bind_rows(gbif_list)
cat("\nTotal GBIF records:", nrow(gbif_raw), "\n")
cat("Records per species:\n")
print(table(gbif_raw$species))

write_csv(gbif_raw, "../data/gbif_raw.csv")
cat("Saved: data/gbif_raw.csv\n")

# =============================================================================
# PART 2: GenBank COI sequences
# =============================================================================
cat("\n=== PART 2: Downloading GenBank sequences ===\n")

fetch_genbank <- function(marker_name, query, retmax = 300) {
    cat("  Marker:", marker_name, "\n")
    cat("  Query:", query, "\n")

    search_res <- entrez_search(db = "nucleotide", term = query, retmax = retmax)
    ids        <- search_res$ids
    cat("  IDs found:", length(ids), "\n")
    if (length(ids) == 0) return(list(fasta = NULL, meta = NULL))

    # Fetch in batches to respect NCBI rate limit (3 req/sec)
    batch_size <- 50
    all_fasta  <- character(0)
    all_meta   <- list()

    for (i in seq(1, length(ids), by = batch_size)) {
        batch_ids <- ids[i:min(i + batch_size - 1, length(ids))]

        fasta_txt <- entrez_fetch(
            db      = "nucleotide",
            id      = batch_ids,
            rettype = "fasta",
            retmode = "text"
        )
        all_fasta <- c(all_fasta, fasta_txt)

        summaries  <- entrez_summary(db = "nucleotide", id = batch_ids)
        batch_meta <- map_dfr(summaries, function(s) {
            tibble(
                accession   = s$accessionversion,
                title       = s$title,
                organism    = s$organism,
                seq_length  = s$slen,
                create_date = s$createdate
            )
        })
        all_meta[[length(all_meta) + 1]] <- batch_meta
        Sys.sleep(0.4)
        cat("  Fetched", min(i + batch_size - 1, length(ids)), "/", length(ids), "\n")
    }

    fasta_path <- paste0("../data/sigaus_", marker_name, "_raw.fasta")
    writeLines(paste(all_fasta, collapse = "\n"), fasta_path)
    cat("  Saved:", fasta_path, "\n")

    list(fasta = fasta_path, meta = bind_rows(all_meta) %>% mutate(marker = marker_name))
}

genbank_results <- imap(GENBANK_QUERIES, function(query, marker) {
    result <- tryCatch(
        fetch_genbank(marker, query),
        error = function(e) { message("ERROR: ", marker, " — ", e$message); NULL }
    )
    Sys.sleep(1)
    result
})

all_meta <- map(genbank_results, "meta") %>%
    keep(~ !is.null(.x)) %>%
    bind_rows()

if (nrow(all_meta) > 0) {
    write_csv(all_meta, "../data/genbank_metadata.csv")
    cat("\nGenBank metadata saved: data/genbank_metadata.csv\n")
    cat("Sequences per marker:\n")
    print(table(all_meta$marker))
}

cat("\n=== Data fetch complete ===\n")
cat("Next: source('02_clean_occurrences.R')\n")
