#!/usr/bin/env Rscript
# Prepare files for METAL meta-analysis
# Add MarkerName column based on configurable format

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
  make_option(c("--markername_format"), type="character", default="chrom:pos:snpid", 
              help="Format for MarkerName column")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Main processing
cat("\n============================================================\n")
cat("Prepare for METAL:", opt$dataset, "\n")
cat("============================================================\n\n")

cat("Input file:", opt$input, "\n")
cat("MarkerName format:", opt$markername_format, "\n\n")

# Read input file
cat("Reading input file...\n")
dt <- fread(opt$input, sep="\t", header=TRUE)
cat("  Loaded", nrow(dt), "variants with", ncol(dt), "columns\n")

# Verify required columns
required_cols <- c("chrom", "pos", "snpid")
missing_cols <- setdiff(required_cols, colnames(dt))

if (length(missing_cols) > 0) {
  stop("ERROR: Missing required columns: ", paste(missing_cols, collapse=", "))
}

# Clean chromosome column - remove "chr" prefix if present
cat("\nCleaning chromosome column...\n")
n_with_chr <- sum(grepl("^chr", dt$chrom, ignore.case=TRUE))
if (n_with_chr > 0) {
  cat("  Removing 'chr' prefix from", n_with_chr, "variants\n")
  dt[, chrom := gsub("^chr", "", chrom, ignore.case=TRUE)]
} else {
  cat("  No 'chr' prefix found\n")
}

# Create MarkerName column based on format
cat("\nCreating MarkerName column...\n")
cat("  Format:", opt$markername_format, "\n")

# Parse format and create markername
if (opt$markername_format == "chrom:pos:snpid") {
  dt[, markername := paste(chrom, pos, snpid, sep=":")]
} else if (opt$markername_format == "chrom:pos") {
  dt[, markername := paste(chrom, pos, sep=":")]
} else if (opt$markername_format == "snpid") {
  dt[, markername := snpid]
} else {
  # Custom format - parse template
  # Replace placeholders with actual column values
  markername_template <- opt$markername_format
  
  # Check for column names in template
  has_chrom <- grepl("chrom", markername_template)
  has_pos <- grepl("pos", markername_template)
  has_snpid <- grepl("snpid", markername_template)
  
  if (has_chrom && has_pos && has_snpid) {
    # Replace each placeholder with column value
    dt[, markername := gsub("chrom", as.character(chrom), markername_template)]
    dt[, markername := gsub("pos", as.character(pos), markername)]
    dt[, markername := gsub("snpid", as.character(snpid), markername)]
  } else {
    # Simple substitution
    markername_template <- gsub("chrom", "%s", markername_template)
    markername_template <- gsub("pos", "%d", markername_template)
    markername_template <- gsub("snpid", "%s", markername_template)
    
    if (has_chrom && has_pos && has_snpid) {
      dt[, markername := sprintf(markername_template, chrom, pos, snpid)]
    } else if (has_chrom && has_pos) {
      dt[, markername := sprintf(markername_template, chrom, pos)]
    } else {
      stop("ERROR: Unsupported markername_format: ", opt$markername_format)
    }
  }
}

cat("  Created MarkerName for", nrow(dt), "variants\n")

# Show examples
cat("\n  MarkerName examples:\n")
example_markers <- head(dt$markername, 5)
for (i in 1:length(example_markers)) {
  cat("    ", example_markers[i], "\n")
}

# Reorder columns - put markername first
current_cols <- colnames(dt)
new_order <- c("markername", setdiff(current_cols, "markername"))
setcolorder(dt, new_order)

cat("\n  Column order:\n")
cat("    ", paste(colnames(dt), collapse=", "), "\n")

# Write output
cat("\nWriting output file...\n")
dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
fwrite(dt, opt$output, sep="\t", quote=FALSE)
cat("  Wrote", nrow(dt), "variants to:", opt$output, "\n")

# Write .done marker
done_content <- paste0(
  "METAL preparation completed for ", opt$dataset, "\n",
  "Input: ", opt$input, "\n",
  "Output: ", opt$output, "\n",
  "Variants: ", nrow(dt), "\n",
  "Columns: ", paste(colnames(dt), collapse=", "), "\n",
  "MarkerName format: ", opt$markername_format, "\n"
)

writeLines(done_content, opt$done)
cat("  Created done marker:", opt$done, "\n")

cat("\n============================================================\n")
cat("METAL preparation completed successfully for", opt$dataset, "\n")
cat("============================================================\n\n")
