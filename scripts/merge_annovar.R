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
  make_option(c("--markername_map"), type="character", help="Markername mapping file (merge_key -> markername)"),
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

# Read markername mapping file (merge_key -> markername)
cat("Reading markername mapping file:", opt$markername_map, "\n")
map_dt <- fread(opt$markername_map, sep="\t", header=FALSE, stringsAsFactors=FALSE)
setnames(map_dt, c("merge_key", "markername"))
cat("Mapping records:", nrow(map_dt), "\n")
cat("Sample mappings (first 5):\n")
print(head(map_dt, 5))
cat("\n")

# Add markername to ANNOVAR annotations by matching Chr:Start:Ref:Alt
cat("Adding markername to ANNOVAR annotations...\n")
annovar_dt[, merge_key := paste(Chr, Start, Ref, Alt, sep=":")]
cat("Sample ANNOVAR merge keys (first 5):", paste(head(annovar_dt$merge_key, 5), collapse=", "), "\n")
annovar_dt <- merge(annovar_dt, map_dt, by="merge_key", all.x=TRUE)
n_matched <- sum(!is.na(annovar_dt$markername))
cat("Matched", n_matched, "out of", nrow(annovar_dt), "ANNOVAR records to markernames\n")
annovar_dt[, merge_key := NULL]  # Remove temporary key
cat("\n")

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

cat("=== ANNOVAR Annotation Merge Complete ===\n")
cat("Total variants:", nrow(annotated_dt), "\n")
cat("Variants with ANNOVAR annotation:", n_with_annotation, "\n")
cat("SNPs without ANNOVAR annotation:", n_without_annotation, "\n\n")

cat("Excel file creation will be handled by separate script\n")

# Close log sinks  
sink(type="message")
sink(type="output")

