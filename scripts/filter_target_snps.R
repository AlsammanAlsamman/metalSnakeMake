#!/usr/bin/env Rscript
# Filter variants based on target SNP list
# If target list not provided or file doesn't exist, copy file unchanged

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--dataset"), type="character", help="Dataset name"),
  make_option(c("--input"), type="character", help="Input TSV file"),
  make_option(c("--output"), type="character", help="Output TSV file"),
  make_option(c("--done"), type="character", help="Done marker file"),
  make_option(c("--log"), type="character", help="Log file"),
  make_option(c("--snp_list"), type="character", default="", help="Target SNP list file (optional)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Main processing
cat("\n============================================================\n")
cat("Filter Target SNPs:", opt$dataset, "\n")
cat("============================================================\n\n")

cat("Input file:", opt$input, "\n")
cat("Target SNP list:", ifelse(opt$snp_list == "", "Not provided", opt$snp_list), "\n\n")

# Read input file
cat("Reading input file...\n")
dt <- fread(opt$input, sep="\t", header=TRUE)
cat("  Loaded", nrow(dt), "variants with", ncol(dt), "columns\n")

initial_count <- nrow(dt)
filtering_applied <- FALSE
target_snps <- NULL

# Check if SNP list is provided and exists
if (opt$snp_list != "" && file.exists(opt$snp_list)) {
  cat("\nReading target SNP list...\n")
  
  # Read SNP list (no header, one column)
  target_snps <- fread(opt$snp_list, header=FALSE, col.names="snpid")
  target_snps[, snpid := trimws(as.character(snpid))]
  target_snps <- target_snps[!is.na(snpid) & snpid != ""]
  
  cat("  Loaded", nrow(target_snps), "target SNPs\n")

  # Empty SNP list means no filtering
  if (nrow(target_snps) == 0) {
    cat("\nTarget SNP list is empty.\n")
    cat("Copying all variants unchanged.\n")
  } else {
  
    # Check if snpid column exists
    if (!"snpid" %in% colnames(dt)) {
      stop("ERROR: 'snpid' column not found in input file. Required for filtering.")
    }
  
    # Filter to target SNPs
    cat("\nFiltering to target SNPs...\n")
    dt_filtered <- dt[snpid %in% target_snps$snpid]
  
    n_kept <- nrow(dt_filtered)
    n_removed <- initial_count - n_kept
  
    cat("  Variants before filtering:", initial_count, "\n")
    cat("  Variants after filtering:", n_kept, "\n")
    cat("  Variants removed:", n_removed, "\n")
    cat("  Retention rate:", round(100 * n_kept / initial_count, 2), "%\n")
  
    # Check for SNPs in target list not found in data
    matched_snps <- sum(target_snps$snpid %in% dt$snpid)
    unmatched_snps <- nrow(target_snps) - matched_snps
  
    cat("\n  Target SNPs found in data:", matched_snps, "/", nrow(target_snps), "\n")
    if (unmatched_snps > 0) {
      cat("  Target SNPs NOT found in data:", unmatched_snps, "\n")
    }
  
    dt <- dt_filtered
    filtering_applied <- TRUE
  }
  
} else {
  if (opt$snp_list != "" && !file.exists(opt$snp_list)) {
    cat("\nWARNING: Target SNP list specified but file not found:", opt$snp_list, "\n")
    cat("Copying all variants unchanged.\n")
  } else {
    cat("\nNo target SNP list provided.\n")
    cat("Copying all variants unchanged.\n")
  }
}

# Write output
cat("\nWriting output file...\n")
dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
fwrite(dt, opt$output, sep="\t", quote=FALSE)
cat("  Wrote", nrow(dt), "variants to:", opt$output, "\n")

# Write .done marker
done_content <- paste0(
  "Filtering completed for ", opt$dataset, "\n",
  "Input: ", opt$input, "\n",
  "Output: ", opt$output, "\n",
  "Input variants: ", initial_count, "\n",
  "Output variants: ", nrow(dt), "\n"
)

if (filtering_applied) {
  done_content <- paste0(
    done_content,
    "Filtering: APPLIED\n",
    "Target SNP list: ", opt$snp_list, "\n",
    "Target SNPs in list: ", nrow(target_snps), "\n",
    "Variants removed: ", initial_count - nrow(dt), "\n",
    "Retention rate: ", round(100 * nrow(dt) / initial_count, 2), "%\n"
  )
} else {
  done_content <- paste0(
    done_content,
    "Filtering: NOT APPLIED (all variants copied)\n"
  )
}

writeLines(done_content, opt$done)
cat("  Created done marker:", opt$done, "\n")

cat("\n============================================================\n")
cat("Filtering completed successfully for", opt$dataset, "\n")
cat("============================================================\n\n")
