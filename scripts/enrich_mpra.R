#!/usr/bin/env Rscript
# Enrich Meta-Analysis with MPRA Functional Data
# Match by chr+pos and aggregate multiple MPRA records per SNP

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(openxlsx)
})

# Parse arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Input meta-analysis enriched file (Step 7 output)"),
  make_option(c("--mpra"), type="character", help="MPRA database file"),
  make_option(c("--output"), type="character", help="Output enriched TSV file"),
  make_option(c("--output_excel"), type="character", help="Output Excel file"),
  make_option(c("--combination"), type="character", help="Combination name"),
  make_option(c("--fdr_threshold"), type="numeric", default=0.05, help="FDR threshold for filtering MPRA [default: 0.05]"),
  make_option(c("--excel_pval"), type="numeric", default=0.005, help="P-value threshold for Excel output [default: 0.005]"),
  make_option(c("--log"), type="character", help="Log file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create log file
log_con <- file(opt$log, open="wt")
sink(log_con, type="output")
sink(log_con, type="message")

cat("=== Enrich Meta-Analysis with MPRA Functional Data ===\n")
cat("Combination:", opt$combination, "\n")
cat("Input meta file:", opt$input, "\n")
cat("MPRA file:", opt$mpra, "\n")
cat("Output TSV:", opt$output, "\n")
cat("Output Excel:", opt$output_excel, "\n")
cat("FDR threshold:", opt$fdr_threshold, "\n")
cat("Excel p-value threshold:", opt$excel_pval, "\n\n")

# Read meta-analysis results
cat("Reading meta-analysis results...\n")
meta_dt <- fread(opt$input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("Meta variants:", nrow(meta_dt), "\n")
cat("Meta columns:", paste(names(meta_dt), collapse=", "), "\n\n")

# Ensure required columns exist
required_cols <- c("chrom", "pos", "p")
missing_cols <- setdiff(required_cols, names(meta_dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns in meta file: ", paste(missing_cols, collapse=", "))
}

# Read MPRA database
cat("Reading MPRA database...\n")
mpra_dt <- fread(opt$mpra, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("MPRA records (raw):", nrow(mpra_dt), "\n")
cat("MPRA columns:", paste(names(mpra_dt), collapse=", "), "\n\n")

# Check MPRA required columns
mpra_required <- c("chr", "pos", "disease", "cellline", "log2FC", "fdr")
mpra_missing <- setdiff(mpra_required, names(mpra_dt))
if (length(mpra_missing) > 0) {
  stop("Missing required columns in MPRA file: ", paste(mpra_missing, collapse=", "))
}

# Filter MPRA for fdr < threshold
cat("Filtering MPRA for fdr <", opt$fdr_threshold, "...\n")
mpra_dt <- mpra_dt[fdr < opt$fdr_threshold]
cat("MPRA records after filtering:", nrow(mpra_dt), "\n\n")

# Convert chromosome column name (chr -> chrom for matching)
setnames(mpra_dt, "chr", "chrom")

# Convert chrom to character for consistent matching
meta_dt[, chrom := as.character(chrom)]
mpra_dt[, chrom := as.character(chrom)]

# Ensure pos is numeric
meta_dt[, pos := as.numeric(pos)]
mpra_dt[, pos := as.numeric(pos)]

# Create matching key
meta_dt[, match_key := paste(chrom, pos, sep=":")]
mpra_dt[, match_key := paste(chrom, pos, sep=":")]

cat("Matching meta SNPs with MPRA by chr:pos...\n")

# Aggregate MPRA records per SNP
cat("Aggregating multiple MPRA records per SNP...\n")
mpra_agg <- mpra_dt[, {
  abs_log2FC <- abs(log2FC[!is.na(log2FC)])
  valid_fdr <- fdr[!is.na(fdr)]
  list(
    mpra_diseases = paste(sort(unique(disease[!is.na(disease) & disease != "None"])), collapse=","),
    mpra_celllines = paste(sort(unique(cellline[!is.na(cellline)])), collapse=","),
    mpra_log2FC_max_abs = if(length(abs_log2FC) > 0) max(abs_log2FC) else NA_real_,
    mpra_log2FC_min_abs = if(length(abs_log2FC) > 0) min(abs_log2FC) else NA_real_,
    mpra_fdr_max = if(length(valid_fdr) > 0) max(valid_fdr) else NA_real_,
    mpra_fdr_min = if(length(valid_fdr) > 0) min(valid_fdr) else NA_real_,
    mpra_n_records = .N
  )
}, by = match_key]

cat("Unique SNPs in MPRA:", nrow(mpra_agg), "\n\n")

# Merge with meta data (left join to keep all meta SNPs)
cat("Merging MPRA data with meta-analysis...\n")
enriched_dt <- merge(meta_dt, mpra_agg, by="match_key", all.x=TRUE)

# Remove temporary match_key column
enriched_dt[, match_key := NULL]

# Count matches
n_with_mpra <- sum(!is.na(enriched_dt$mpra_n_records))
n_without_mpra <- sum(is.na(enriched_dt$mpra_n_records))

cat("SNPs with MPRA data:", n_with_mpra, "\n")
cat("SNPs without MPRA data:", n_without_mpra, "\n\n")

# Write full TSV output
cat("Writing TSV output...\n")
fwrite(enriched_dt, opt$output, sep="\t", quote=FALSE, na="NA")
cat("TSV output written:", opt$output, "\n\n")

# Filter for Excel output using threshold from config
cat("Filtering variants for Excel output (p <", opt$excel_pval, ")...\n")
excel_dt <- enriched_dt[p < opt$excel_pval]
cat("Variants for Excel:", nrow(excel_dt), "\n\n")

# Create Excel output with formatting
cat("Creating Excel output with formatting...\n")
wb <- createWorkbook()
addWorksheet(wb, "MPRA_Enriched")

# Write data
writeData(wb, "MPRA_Enriched", excel_dt)

# Format header
headerStyle <- createStyle(
  fontSize = 12,
  fontColour = "#FFFFFF",
  halign = "center",
  fgFill = "#4F81BD",
  border = "TopBottom",
  borderColour = "#4F81BD",
  textDecoration = "bold"
)
addStyle(wb, sheet = "MPRA_Enriched", headerStyle, rows = 1, cols = 1:ncol(excel_dt), gridExpand = TRUE)

# Freeze first row
freezePane(wb, "MPRA_Enriched", firstRow = TRUE)

# Highlight SNPs with MPRA data
if (n_with_mpra > 0) {
  mpra_rows <- which(!is.na(excel_dt$mpra_n_records)) + 1  # +1 for header
  if (length(mpra_rows) > 0) {
    mpra_style <- createStyle(fgFill = "#E6F3FF")
    addStyle(wb, sheet = "MPRA_Enriched", mpra_style, rows = mpra_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE)
    cat("Highlighted", length(mpra_rows), "rows with MPRA data\n")
  }
}

# Highlight genome-wide significant SNPs (p < 5e-8) in different color
sig_rows <- which(excel_dt$p < 5e-8) + 1  # +1 for header
if (length(sig_rows) > 0) {
  sig_style <- createStyle(fgFill = "#FFE6E6", textDecoration = "bold")
  addStyle(wb, sheet = "MPRA_Enriched", sig_style, rows = sig_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE, stack = TRUE)
  cat("Highlighted", length(sig_rows), "rows with p < 5e-8\n")
}

# Auto-size columns
setColWidths(wb, sheet = "MPRA_Enriched", cols = 1:ncol(excel_dt), widths = "auto")

# Save Excel file
saveWorkbook(wb, opt$output_excel, overwrite = TRUE)
cat("Excel output written:", opt$output_excel, "\n\n")

cat("=== MPRA Enrichment Complete ===\n")
cat("Total variants:", nrow(enriched_dt), "\n")
cat("Variants with MPRA data:", n_with_mpra, "\n")
cat("Variants in Excel (p < 5e-3):", nrow(excel_dt), "\n")

# Close log
sink(type="message")
sink(type="output")
close(log_con)
