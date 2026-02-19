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
  make_option(c("--studies"), type="character", default="", help="Comma-separated list of studies used in meta-analysis"),
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
cat("Studies:", opt$studies, "\n\n")

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

# Add meta-level OR if missing (OR = exp(beta))
if (!"or" %in% names(annotated_dt) && "beta" %in% names(annotated_dt)) {
  cat("Adding meta OR from beta (or = exp(beta))...\n")
  annotated_dt[, or := suppressWarnings(exp(as.numeric(beta)))]
} else if ("or" %in% names(annotated_dt) && "beta" %in% names(annotated_dt)) {
  cat("Filling missing meta OR values from beta...\n")
  annotated_dt[is.na(or) & !is.na(beta), or := suppressWarnings(exp(as.numeric(beta)))]
}

# Merge individual study-level data used in meta-analysis
study_names <- trimws(unlist(strsplit(opt$studies, ",", fixed=TRUE)))
study_names <- study_names[nchar(study_names) > 0]

if (length(study_names) > 0) {
  cat("Merging individual study data for:", paste(study_names, collapse=", "), "\n")

  for (study in study_names) {
    study_file <- file.path("results", "04_metal_ready", paste0(study, ".tsv"))

    if (!file.exists(study_file)) {
      cat("  WARNING: Study file not found:", study_file, "\n")
      next
    }

    cat("  Reading", study, "from", study_file, "...\n")
    study_dt <- fread(study_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

    # Ensure markername exists for joining
    if (!"markername" %in% names(study_dt)) {
      if (all(c("chrom", "pos", "snpid") %in% names(study_dt))) {
        study_dt[, markername := paste(chrom, pos, snpid, sep=":")]
      } else {
        cat("  WARNING: Cannot create markername for", study, "(missing markername/chrom+pos+snpid), skipping\n")
        next
      }
    }

    # Add study OR if missing
    if (!"or" %in% names(study_dt) && "beta" %in% names(study_dt)) {
      study_dt[, or := suppressWarnings(exp(as.numeric(beta)))]
    } else if ("or" %in% names(study_dt) && "beta" %in% names(study_dt)) {
      study_dt[is.na(or) & !is.na(beta), or := suppressWarnings(exp(as.numeric(beta)))]
    }

    # Keep key study columns only
    keep_cols <- intersect(c("markername", "beta", "se", "p", "eaf", "or"), names(study_dt))
    if (!"markername" %in% keep_cols) {
      cat("  WARNING: markername still missing in", study, "after processing, skipping\n")
      next
    }
    study_subset <- study_dt[, ..keep_cols]

    # Prefix study columns
    rename_cols <- setdiff(names(study_subset), "markername")
    if (length(rename_cols) > 0) {
      setnames(study_subset, rename_cols, paste0(study, "_", rename_cols))
    }

    # Merge into annotated table
    n_before <- ncol(annotated_dt)
    annotated_dt <- merge(annotated_dt, study_subset, by="markername", all.x=TRUE)
    n_added <- ncol(annotated_dt) - n_before
    cat("  Added", n_added, "columns for", study, "\n")
  }

  cat("Individual study data merge complete.\n\n")
}

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

