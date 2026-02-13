#!/usr/bin/env Rscript
# Enrich Meta-Analysis with Study-Specific EAF
# Align EAF values to meta effect allele and calculate meta EAF

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(openxlsx)
})

# Parse arguments
option_list <- list(
  make_option(c("--meta"), type="character", help="Input meta-analysis file (Step 6 standardized output)"),
  make_option(c("--output"), type="character", help="Output enriched TSV file"),
  make_option(c("--output_excel"), type="character", help="Output enriched Excel file"),
  make_option(c("--combination"), type="character", help="Study combination name"),
  make_option(c("--studies"), type="character", help="Comma-separated list of studies"),
  make_option(c("--log"), type="character", help="Log file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create log file
log_con <- file(opt$log, open="wt")
sink(log_con, type="output")
sink(log_con, type="message")

cat("=== Enrich Meta-Analysis with Study-Specific EAF ===\n")
cat("Combination:", opt$combination, "\n")
cat("Studies:", opt$studies, "\n")
cat("Input meta file:", opt$meta, "\n")
cat("Output TSV:", opt$output, "\n")
cat("Output Excel:", opt$output_excel, "\n\n")

# Parse studies
studies <- unlist(strsplit(opt$studies, ","))
cat("Number of studies:", length(studies), "\n")
cat("Study list:", paste(studies, collapse=", "), "\n\n")

# Read meta-analysis results
cat("Reading meta-analysis results...\n")
meta_dt <- fread(opt$meta, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("Meta variants:", nrow(meta_dt), "\n")
cat("Meta columns:", paste(names(meta_dt), collapse=", "), "\n\n")

# Ensure required columns exist
required_cols <- c("markername", "ea", "nea")
missing_cols <- setdiff(required_cols, names(meta_dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns in meta file: ", paste(missing_cols, collapse=", "))
}

# Standardize allele columns to uppercase for matching
meta_dt[, ea := toupper(ea)]
meta_dt[, nea := toupper(nea)]

# Vectorized function to align EAF to meta effect allele
align_eaf_vectorized <- function(dt) {
  # Initialize aligned_eaf as NA
  dt[, aligned_eaf := NA_real_]
  
  # Case 1: No flip needed - study EA matches meta EA
  dt[!is.na(study_ea) & !is.na(ea) & study_ea == ea, aligned_eaf := study_eaf]
  n_case1 <- dt[!is.na(study_ea) & !is.na(ea) & study_ea == ea, .N]
  
  # Case 2: Flip needed - study EA matches meta NEA (study's effect allele is meta's non-effect allele)
  dt[!is.na(study_ea) & !is.na(nea) & study_ea == nea & is.na(aligned_eaf), aligned_eaf := 1 - study_eaf]
  n_case2 <- dt[!is.na(study_ea) & !is.na(nea) & study_ea == nea, .N] - n_case1
  
  # Case 3: Flip needed - study NEA matches meta EA (study's non-effect allele is meta's effect allele)
  dt[!is.na(study_oa) & !is.na(ea) & study_oa == ea & is.na(aligned_eaf), aligned_eaf := 1 - study_eaf]
  n_case3 <- dt[!is.na(study_oa) & !is.na(ea) & study_oa == ea & is.na(aligned_eaf), .N]
  
  return(list(
    aligned_eaf = dt$aligned_eaf,
    n_case1 = n_case1,
    n_case2 = max(0, n_case2),
    n_case3 = n_case3
  ))
}

# Read and merge study-specific EAF data
for (study in studies) {
  cat("Processing study:", study, "\n")
  
  # Read study data from Step 4 (metal_ready)
  study_file <- file.path("results", "04_metal_ready", paste0(study, ".tsv"))
  
  if (!file.exists(study_file)) {
    cat("  WARNING: Study file not found:", study_file, "\n")
    cat("  Setting all EAF values to NA for", study, "\n")
    meta_dt[, paste0("eaf_", study) := NA_real_]
    next
  }
  
  study_dt <- fread(study_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  cat("  Study variants:", nrow(study_dt), "\n")
  
  # Check required columns
  study_required <- c("markername", "ea", "nea", "eaf")
  study_missing <- setdiff(study_required, names(study_dt))
  if (length(study_missing) > 0) {
    cat("  WARNING: Missing columns in study file:", paste(study_missing, collapse=", "), "\n")
    meta_dt[, paste0("eaf_", study) := NA_real_]
    next
  }
  
  # Standardize allele columns
  study_dt[, ea := toupper(ea)]
  study_dt[, nea := toupper(nea)]
  
  # Remove duplicate markernames (keep first occurrence)
  n_before_dedup <- nrow(study_dt)
  study_dt <- study_dt[!duplicated(markername)]
  n_after_dedup <- nrow(study_dt)
  if (n_before_dedup > n_after_dedup) {
    cat("  Removed", n_before_dedup - n_after_dedup, "duplicate markernames\n")
  }
  
  # Select and rename columns
  study_dt <- study_dt[, .(markername, 
                           study_ea = ea, 
                           study_oa = nea, 
                           study_eaf = eaf)]
  
  # Merge with meta data
  merged <- merge(meta_dt[, .(markername, ea, nea)], 
                  study_dt, 
                  by = "markername", 
                  all.x = TRUE)
  
  # Align EAF to meta effect allele (vectorized)
  align_result <- align_eaf_vectorized(merged)
  
  # Add aligned EAF to the merged data
  merged[, aligned_eaf := align_result$aligned_eaf]
  
  # Merge aligned EAF back to meta_dt by markername (preserves row order)
  eaf_col <- paste0("eaf_", study)
  meta_dt[merged, (eaf_col) := i.aligned_eaf, on = "markername"]
  
  # Report statistics
  n_aligned <- sum(!is.na(align_result$aligned_eaf))
  n_flipped <- align_result$n_case2 + align_result$n_case3
  cat("  Aligned variants:", n_aligned, "\n")
  cat("  No flip (EA match):", align_result$n_case1, "\n")
  cat("  Flipped EAF (1-EAF):", n_flipped, "\n")
}

cat("\n")

# Calculate meta EAF as mean of aligned study EAFs (excluding NAs)
cat("Calculating meta EAF as mean of aligned study EAFs...\n")
eaf_cols <- paste0("eaf_", studies)

# Create matrix of study EAFs
eaf_matrix <- as.matrix(meta_dt[, ..eaf_cols])

# Calculate mean excluding NAs and count of non-NA studies
meta_dt[, eaf := rowMeans(eaf_matrix, na.rm = TRUE)]
meta_dt[, n_studies_with_eaf := rowSums(!is.na(eaf_matrix))]

# Set eaf to NA if no studies have data
meta_dt[n_studies_with_eaf == 0, eaf := NA_real_]

cat("Variants with EAF from all studies:", sum(meta_dt$n_studies_with_eaf == length(studies)), "\n")
cat("Variants with EAF from some studies:", sum(meta_dt$n_studies_with_eaf > 0 & meta_dt$n_studies_with_eaf < length(studies)), "\n")
cat("Variants with no EAF data:", sum(meta_dt$n_studies_with_eaf == 0), "\n\n")

# Reorder columns: meta columns first, then eaf, then study-specific EAFs, then n_studies_with_eaf
meta_cols <- setdiff(names(meta_dt), c(eaf_cols, "eaf", "n_studies_with_eaf"))
setcolorder(meta_dt, c(meta_cols, "eaf", eaf_cols, "n_studies_with_eaf"))

# Write TSV output
cat("Writing TSV output...\n")
fwrite(meta_dt, opt$output, sep="\t", quote=FALSE, na="NA")
cat("TSV output written:", opt$output, "\n\n")

# Filter for Excel output (p < 5e-5)
cat("Filtering variants for Excel output (p < 5e-5)...\n")
excel_dt <- meta_dt[p < 5e-5]
cat("Variants for Excel:", nrow(excel_dt), "\n\n")

# Create Excel output with formatting
cat("Creating Excel output with formatting...\n")

wb <- createWorkbook()
addWorksheet(wb, "Meta_Analysis")

# Write data
writeData(wb, sheet=1, excel_dt, startRow=1, startCol=1, headerStyle=createStyle(textDecoration="bold"))

# Apply formatting
# Scientific notation for p-value column
if ("p" %in% names(excel_dt)) {
  p_col <- which(names(excel_dt) == "p")
  style_scientific <- createStyle(numFmt="0.00E+00")
  addStyle(wb, sheet=1, style=style_scientific, rows=2:(nrow(excel_dt)+1), cols=p_col, gridExpand=TRUE)
}

# Comma format for position
if ("pos" %in% names(excel_dt)) {
  pos_col <- which(names(excel_dt) == "pos")
  style_comma <- createStyle(numFmt="#,##0")
  addStyle(wb, sheet=1, style=style_comma, rows=2:(nrow(excel_dt)+1), cols=pos_col, gridExpand=TRUE)
}

# 2 decimal places for numeric columns (including EAF columns)
numeric_cols <- c("beta", "se", "or", "eaf", eaf_cols)
numeric_cols <- intersect(numeric_cols, names(excel_dt))
style_decimal <- createStyle(numFmt="0.00")
for (col_name in numeric_cols) {
  col_idx <- which(names(excel_dt) == col_name)
  addStyle(wb, sheet=1, style=style_decimal, rows=2:(nrow(excel_dt)+1), cols=col_idx, gridExpand=TRUE)
}

# Highlight p < 5e-8 (if p column exists)
if ("p" %in% names(excel_dt)) {
  p_col <- which(names(excel_dt) == "p")
  sig_rows <- which(excel_dt$p < 5e-8) + 1  # +1 for header row
  if (length(sig_rows) > 0) {
    style_highlight <- createStyle(fgFill="#FFB3B3")
    for (row_idx in sig_rows) {
      addStyle(wb, sheet=1, style=style_highlight, rows=row_idx, cols=1:ncol(excel_dt), gridExpand=TRUE)
    }
    cat("Highlighted", length(sig_rows), "rows with p < 5e-8\n")
  }
}

# Auto-size columns
setColWidths(wb, sheet=1, cols=1:ncol(excel_dt), widths="auto")

# Save Excel file
saveWorkbook(wb, opt$output_excel, overwrite=TRUE)
cat("Excel output written:", opt$output_excel, "\n\n")

cat("=== EAF Enrichment Complete ===\n")
cat("Total variants:", nrow(meta_dt), "\n")
cat("Study-specific EAF columns added:", length(eaf_cols), "\n")
cat("Meta EAF column name: eaf\n")

# Close log
sink(type="message")
sink(type="output")
close(log_con)
