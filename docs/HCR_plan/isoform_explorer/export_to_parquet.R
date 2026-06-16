#!/usr/bin/env Rscript
# export_to_parquet.R
#
# One-time R export: reads .qs constitutive-sequence_list files for all three
# species and writes intermediate CSVs that the companion Python script
# (convert_to_parquet.py) then converts to Parquet + zstd.
#
# Run from the repo root:
#   Rscript docs/HCR_plan/isoform_explorer/export_to_parquet.R
#
# After verifying the Parquet files are correct, clean up:
#   rm /tmp/isoform_explorer_csv/*.csv
#   rm data/isoform_explorer/{Hsap,Mmus,Dmel}/*.qs
#   rm data/isoform_explorer/{Hsap,Mmus,Dmel}/*_summary.csv

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(purrr)
})

DATA_DIR <- "data/isoform_explorer"
OUT_DIR  <- "/tmp/isoform_explorer_csv"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── TSL cutoff string from min TSL value (mammals only) ──────────────────────
derive_tsl_cutoff <- function(min_tsl) {
  if (is.na(min_tsl) || min_tsl > 5) {
    "all (no TSL)"
  } else if (min_tsl <= 2) {
    "TSL<=2"
  } else {
    paste0("TSL==", as.integer(min_tsl))
  }
}

# ── Per-species export ────────────────────────────────────────────────────────
export_species <- function(sp, qs_path, is_mammal) {
  cat("\n=== Processing", sp, "===\n")
  cat("Reading:", qs_path, "\n")
  lst <- qread(qs_path)
  cat("Loaded", length(lst), "genes.\n")

  # Step 0 — verify transcript_support_level is present for mammals
  if (is_mammal) {
    sample_exon_cols <- names(lst[[1]][["exon_intervals"]])
    if (!"transcript_support_level" %in% sample_exon_cols) {
      stop(
        sp, ": transcript_support_level column MISSING from exon_intervals!\n",
        "  Available columns: ", paste(sample_exon_cols, collapse = ", "), "\n",
        "  Cannot build TSL-coloured plot. Re-join TSL from the GTF before exporting."
      )
    }
    tsl_sample <- unique(lst[[1]][["exon_intervals"]][["transcript_support_level"]])
    cat("  [OK] transcript_support_level present; sample values:", paste(tsl_sample, collapse = ", "), "\n")
  }

  # ── Summary tibble (one row per gene) ──────────────────────────────────────
  cat("Building summary...\n")
  summary_df <- map_dfr(lst, function(g) {
    tsl_cutoff   <- if (is_mammal) derive_tsl_cutoff(g$min_transcript_support_level) else "all isoforms (no TSL)"
    isoform_count <- if (is_mammal) g$tsl_isoform_count else g$isoform_count
    seq_val <- g$constitutive_seq_Ncollapsed
    tibble(
      gene_id                   = g$gene_id,
      gene_name                 = g$gene_name,
      gene_biotype              = g$gene_biotype,
      strand                    = g$strand,
      isoform_count             = as.integer(isoform_count),
      constitutive_exist        = as.logical(g$constitutive_exist),
      constitutive_length       = as.integer(g$constitutive_length),
      constitutive_seq_Ncollapsed = if (is.null(seq_val) || length(seq_val) == 0) NA_character_ else as.character(seq_val),
      tsl_cutoff                = tsl_cutoff
    )
  }, .progress = TRUE)

  summary_path <- file.path(OUT_DIR, paste0(sp, "_summary.csv"))
  write.csv(summary_df, summary_path, row.names = FALSE, na = "")
  cat("  Summary:", nrow(summary_df), "rows ->", summary_path, "\n")

  # ── Intervals tibble (long format, drives the plot) ────────────────────────
  cat("Building intervals...\n")
  intervals_df <- map_dfr(lst, function(g) {
    gid    <- g$gene_id
    strand <- g$strand

    # Exon rows
    ei <- g$exon_intervals
    if (is_mammal) {
      tsl_col     <- as.integer(ei$transcript_support_level)
      tsl_included <- ei$transcript_id %in% g$tsl_filtered_transcript_ids
    } else {
      tsl_col      <- rep(NA_integer_, nrow(ei))
      tsl_included <- rep(TRUE, nrow(ei))
    }

    exon_rows <- tibble(
      gene_id                  = gid,
      transcript_id            = ei$transcript_id,
      feature                  = "exon",
      start                    = as.integer(ei$start),
      end                      = as.integer(ei$end),
      strand                   = strand,
      transcript_support_level = tsl_col,
      tsl_included             = tsl_included
    )

    # Constitutive rows (may be 0 rows if no constitutive seq)
    ci <- g$constitutive_intervals
    if (!is.null(ci) && nrow(ci) > 0) {
      const_rows <- tibble(
        gene_id                  = gid,
        transcript_id            = "CONSTITUTIVE",
        feature                  = "constitutive",
        start                    = as.integer(ci$start),
        end                      = as.integer(ci$end),
        strand                   = strand,
        transcript_support_level = NA_integer_,
        tsl_included             = TRUE
      )
      bind_rows(exon_rows, const_rows)
    } else {
      exon_rows
    }
  }, .progress = TRUE)

  # Sort by gene_id so pyarrow predicate-pushdown can skip row groups
  intervals_df <- intervals_df %>% arrange(gene_id)

  intervals_path <- file.path(OUT_DIR, paste0(sp, "_intervals.csv"))
  write.csv(intervals_df, intervals_path, row.names = FALSE, na = "")
  cat("  Intervals:", nrow(intervals_df), "rows ->", intervals_path, "\n")
}

# ── Run all three species ─────────────────────────────────────────────────────
export_species(
  "Hsap",
  file.path(DATA_DIR, "Hsap/Hsap_ensembl_v113_tsl_proteincoding_noncoding_constitutive-sequence_list.qs"),
  is_mammal = TRUE
)

export_species(
  "Mmus",
  file.path(DATA_DIR, "Mmus/Mmus_ensembl_v113_tsl_proteincoding_noncoding_constitutive-sequence_list.qs"),
  is_mammal = TRUE
)

export_species(
  "Dmel",
  file.path(DATA_DIR, "Dmel/Dmel_ensembl_v113_proteincoding_noncoding_constitutive-sequence_list.qs"),
  is_mammal = FALSE
)

cat("\n=== Done. Intermediate CSVs in", OUT_DIR, "===\n")
cat("Next step:\n")
cat("  /Users/jefflee/miniforge3/envs/probedesign/bin/python",
    "docs/HCR_plan/isoform_explorer/convert_to_parquet.py\n")
