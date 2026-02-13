#!/usr/bin/env Rscript
# Merge ANNOVAR Annotations with Enriched Meta-Analysis File

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(openxlsx)
})

# Parse arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Input meta-analysis enriched file"),
  make_option(c("--annovar"), type="character", help="ANNOVAR annotation output file"),
  make_option(c("--output"), type="character", help="Output annotated TSV file"),
  make_option(c("--output_excel"), type="character", help="Output Excel file"),
  make_option(c("--excel_pval"), type="numeric", help="P-value threshold for Excel output"),
  make_option(c("--combination"), type="character", help="Combination name"),
  make_option(c("--log"), type="character", help="Log file path")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Append to existing log
log_con <- file(opt$log, open="at")
sink(log_con, type="output", append=TRUE)
sink(log_con, type="message", append=TRUE)

cat("\n=== Merging ANNOVAR Annotations ===\n")
cat("Input meta file:", opt$input, "\n")
cat("ANNOVAR file:", opt$annovar, "\n")
cat("Output TSV:", opt$output, "\n")
cat("Output Excel:", opt$output_excel, "\n")
cat("Excel p-value threshold:", opt$excel_pval, "\n\n")

# Read enriched meta-analysis file
cat("Reading enriched meta-analysis file...\n")
meta_dt <- fread(opt$input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("Meta variants:", nrow(meta_dt), "\n")
cat("Meta columns:", paste(names(meta_dt), collapse=", "), "\n\n")

# Read ANNOVAR annotations
cat("Reading ANNOVAR annotations...\n")
annovar_dt <- fread(opt$annovar, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("ANNOVAR records:", nrow(annovar_dt), "\n")
cat("ANNOVAR columns:", paste(names(annovar_dt), collapse=", "), "\n\n")

# ANNOVAR output includes the key in "Otherinfo1" column (last column we added)
# Rename for clarity
if ("Otherinfo1" %in% names(annovar_dt)) {
  setnames(annovar_dt, "Otherinfo1", "markername")
}

# Check if markername exists
if (!"markername" %in% names(annovar_dt)) {
  stop("ERROR: markername column not found in ANNOVAR output")
}

# Select ANNOVAR annotation columns (exclude input columns: Chr, Start, End, Ref, Alt)
annovar_cols <- setdiff(names(annovar_dt), c("Chr", "Start", "End", "Ref", "Alt"))
annovar_dt <- annovar_dt[, ..annovar_cols]

cat("ANNOVAR annotation columns to add:", paste(annovar_cols[annovar_cols != "markername"], collapse=", "), "\n\n")

# Merge ANNOVAR annotations with meta data by markername
cat("Merging ANNOVAR annotations with meta-analysis...\n")
annotated_dt <- merge(meta_dt, annovar_dt, by="markername", all.x=TRUE)

# Count matches
n_with_annotation <- sum(!is.na(annotated_dt$Func.refGene))
n_without_annotation <- sum(is.na(annotated_dt$Func.refGene))

cat("SNPs with ANNOVAR annotation:", n_with_annotation, "\n")
cat("SNPs without ANNOVAR annotation:", n_without_annotation, "\n\n")

# Write full TSV output
cat("Writing TSV output...\n")
fwrite(annotated_dt, opt$output, sep="\t", quote=FALSE, na="NA")
cat("TSV output written:", opt$output, "\n\n")

# Filter for Excel output using threshold from config
cat("Filtering variants for Excel output (p <", opt$excel_pval, ")...\n")
excel_dt <- annotated_dt[p < opt$excel_pval]
cat("Variants for Excel:", nrow(excel_dt), "\n\n")

# Create Excel output with formatting
cat("Creating Excel output with formatting...\n")
wb <- createWorkbook()
addWorksheet(wb, "ANNOVAR_Annotated")

# Write data
writeData(wb, "ANNOVAR_Annotated", excel_dt)

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
addStyle(wb, sheet = "ANNOVAR_Annotated", headerStyle, rows = 1, cols = 1:ncol(excel_dt), gridExpand = TRUE)

# Freeze first row
freezePane(wb, "ANNOVAR_Annotated", firstRow = TRUE)

# Highlight SNPs with ANNOVAR annotations
if (n_with_annotation > 0 && nrow(excel_dt) > 0) {
  annotated_rows <- which(!is.na(excel_dt$Func.refGene)) + 1  # +1 for header
  if (length(annotated_rows) > 0) {
    annotated_style <- createStyle(fgFill = "#E6F3FF")
    addStyle(wb, sheet = "ANNOVAR_Annotated", annotated_style, rows = annotated_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE)
    cat("Highlighted", length(annotated_rows), "rows with ANNOVAR annotation\n")
  }
}

# Highlight genome-wide significant SNPs (p < 5e-8) in different color
sig_rows <- which(excel_dt$p < 5e-8) + 1  # +1 for header
if (length(sig_rows) > 0) {
  sig_style <- createStyle(fgFill = "#FFE6E6", textDecoration = "bold")
  addStyle(wb, sheet = "ANNOVAR_Annotated", sig_style, rows = sig_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE, stack = TRUE)
  cat("Highlighted", length(sig_rows), "rows with p < 5e-8\n")
}

# Auto-size columns
setColWidths(wb, sheet = "ANNOVAR_Annotated", cols = 1:ncol(excel_dt), widths = "auto")

# Save Excel file
saveWorkbook(wb, opt$output_excel, overwrite = TRUE)
cat("Excel output written:", opt$output_excel, "\n\n")

cat("=== ANNOVAR Annotation Merge Complete ===\n")
cat("Total variants:", nrow(annotated_dt), "\n")
cat("Variants with ANNOVAR annotation:", n_with_annotation, "\n")
cat("Variants in Excel (p <", opt$excel_pval, "):", nrow(excel_dt), "\n")

# Close log
sink(type="message")
sink(type="output")
close(log_con)
