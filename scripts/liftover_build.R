#!/usr/bin/env Rscript
# LiftOver genome coordinates to target build
# If dataset already in target build, copy unchanged
# If different build, perform liftOver and update chrom/pos

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
  make_option(c("--source_build"), type="character", help="Source genome build (e.g., hg19)"),
  make_option(c("--target_build"), type="character", help="Target genome build (e.g., hg38)"),
  make_option(c("--chain_file"), type="character", default="NA", help="LiftOver chain file path (NA if not needed)"),
  make_option(c("--liftover_bin"), type="character", help="Path to liftOver binary")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Main processing
cat("\n============================================================\n")
cat("LiftOver Genome Coordinates:", opt$dataset, "\n")
cat("============================================================\n\n")

cat("Input file:", opt$input, "\n")
cat("Source build:", opt$source_build, "\n")
cat("Target build:", opt$target_build, "\n\n")

# Read input file
cat("Reading input file...\n")
dt <- fread(opt$input, sep="\t", header=TRUE)
cat("  Loaded", nrow(dt), "variants with", ncol(dt), "columns\n")

initial_count <- nrow(dt)
liftover_applied <- FALSE
lifted_count <- 0
unmapped_count <- 0

# Check if liftOver is needed
if (opt$source_build == opt$target_build) {
  cat("\nSource and target builds are identical (", opt$target_build, ")\n")
  cat("Copying all variants unchanged.\n")
  
} else {
  cat("\nBuild conversion required:", opt$source_build, "->", opt$target_build, "\n")
  
  # Verify required columns
  if (!all(c("chrom", "pos") %in% colnames(dt))) {
    stop("ERROR: 'chrom' and 'pos' columns required for liftOver")
  }
  
  # Check chain file
  if (opt$chain_file == "" || opt$chain_file == "NA" || !file.exists(opt$chain_file)) {
    stop("ERROR: Chain file not found: ", opt$chain_file)
  }
  cat("Chain file:", opt$chain_file, "\n")
  
  # Check liftOver binary - look in PATH if not absolute path
  if (file.exists(opt$liftover_bin)) {
    liftover_path <- opt$liftover_bin
  } else {
    # Try to find binary in PATH using Sys.which
    liftover_path <- Sys.which(opt$liftover_bin)
    if (liftover_path == "") {
      stop("ERROR: liftOver binary not found in PATH: ", opt$liftover_bin, 
           "\n  Make sure the module is loaded or provide full path")
    }
  }
  cat("LiftOver binary:", liftover_path, "\n\n")
  
  # Validate chromosome format for source build
  cat("Validating chromosome format for source build...\n")
  unique_chroms <- unique(dt$chrom)
  cat("  Unique chromosomes:", paste(head(unique_chroms, 10), collapse=", "), 
      ifelse(length(unique_chroms) > 10, "...", ""), "\n")
  
  # Standardize chromosome names - add "chr" prefix if missing
  has_chr_prefix <- any(grepl("^chr", dt$chrom))
  if (!has_chr_prefix) {
    cat("  Adding 'chr' prefix to chromosome names\n")
    dt[, chrom := paste0("chr", chrom)]
  }
  
  # Create temporary BED file (0-based for liftOver)
  temp_dir <- tempdir()
  input_bed <- file.path(temp_dir, paste0(opt$dataset, "_input.bed"))
  output_bed <- file.path(temp_dir, paste0(opt$dataset, "_lifted.bed"))
  unmapped_bed <- file.path(temp_dir, paste0(opt$dataset, "_unmapped.bed"))
  
  cat("\nPreparing BED file for liftOver...\n")
  # BED format: chrom, start (0-based), end (1-based), name
  # For SNPs: start = pos - 1, end = pos
  bed_df <- data.frame(
    chrom = dt$chrom,
    start = as.integer(dt$pos - 1),  # Convert to 0-based, ensure integer
    end = as.integer(dt$pos),        # Ensure integer
    name = paste0("var_", 1:nrow(dt)),  # Unique identifier for each variant
    stringsAsFactors = FALSE
  )
  
  # Write BED file with scipen to avoid scientific notation
  options(scipen = 999)  # Prevent scientific notation
  fwrite(bed_df, input_bed, sep="\t", col.names=FALSE, quote=FALSE, scipen=999)
  options(scipen = 0)  # Reset
  cat("  Created BED file with", nrow(bed_df), "variants\n")
  
  # Run liftOver
  cat("\nRunning liftOver...\n")
  liftover_cmd <- sprintf(
    "%s %s %s %s %s",
    liftover_path,
    input_bed,
    opt$chain_file,
    output_bed,
    unmapped_bed
  )
  
  system_result <- system(liftover_cmd, intern=FALSE)
  
  if (system_result != 0) {
    stop("ERROR: liftOver command failed with exit code ", system_result)
  }
  
  # Read liftOver results
  if (file.exists(output_bed) && file.info(output_bed)$size > 0) {
    lifted_bed <- fread(output_bed, sep="\t", header=FALSE, 
                        col.names=c("chrom", "start", "end", "name"))
    lifted_count <- nrow(lifted_bed)
    cat("  Successfully lifted", lifted_count, "variants\n")
    
    # Convert back to 1-based position
    lifted_bed[, pos := end]
    
    # Extract original row index from name
    lifted_bed[, row_idx := as.integer(sub("var_", "", name))]
    
    # Merge lifted coordinates back to original data
    dt[, row_idx := 1:.N]
    dt_lifted <- merge(dt, lifted_bed[, .(row_idx, chrom_new=chrom, pos_new=pos)], 
                       by="row_idx", all.x=FALSE)  # Keep only lifted variants
    
    # Update chrom and pos
    dt_lifted[, chrom := chrom_new]
    dt_lifted[, pos := pos_new]
    dt_lifted[, c("row_idx", "chrom_new", "pos_new") := NULL]
    
    dt <- dt_lifted
    
  } else {
    stop("ERROR: liftOver produced no output")
  }
  
  # Check for unmapped variants
  if (file.exists(unmapped_bed) && file.info(unmapped_bed)$size > 0) {
    unmapped_lines <- readLines(unmapped_bed)
    # Count actual variant lines (skip comment lines starting with #)
    unmapped_count <- sum(!grepl("^#", unmapped_lines))
    cat("  Unmapped variants (removed):", unmapped_count, "\n")
  }
  
  cat("\n  Variants before liftOver:", initial_count, "\n")
  cat("  Variants after liftOver:", nrow(dt), "\n")
  cat("  Variants removed:", initial_count - nrow(dt), "\n")
  cat("  Success rate:", round(100 * nrow(dt) / initial_count, 2), "%\n")
  
  # Validate chromosome format for target build
  cat("\nValidating chromosome format for target build...\n")
  unique_chroms_new <- unique(dt$chrom)
  cat("  Unique chromosomes after liftOver:", 
      paste(head(unique_chroms_new, 10), collapse=", "),
      ifelse(length(unique_chroms_new) > 10, "...", ""), "\n")
  
  # Check for non-standard chromosomes
  standard_chroms <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")
  non_standard <- setdiff(unique_chroms_new, standard_chroms)
  if (length(non_standard) > 0) {
    cat("  WARNING: Non-standard chromosomes detected:", 
        paste(head(non_standard, 5), collapse=", "), "\n")
  }
  
  # Clean up temporary files
  unlink(c(input_bed, output_bed, unmapped_bed))
  
  liftover_applied <- TRUE
}

# Write output
cat("\nWriting output file...\n")
dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
fwrite(dt, opt$output, sep="\t", quote=FALSE)
cat("  Wrote", nrow(dt), "variants to:", opt$output, "\n")

# Write .done marker
done_content <- paste0(
  "LiftOver completed for ", opt$dataset, "\n",
  "Input: ", opt$input, "\n",
  "Output: ", opt$output, "\n",
  "Source build: ", opt$source_build, "\n",
  "Target build: ", opt$target_build, "\n",
  "Input variants: ", initial_count, "\n",
  "Output variants: ", nrow(dt), "\n"
)

if (liftover_applied) {
  done_content <- paste0(
    done_content,
    "LiftOver: APPLIED\n",
    "Chain file: ", opt$chain_file, "\n",
    "Variants lifted: ", lifted_count, "\n",
    "Variants unmapped: ", unmapped_count, "\n",
    "Variants removed: ", initial_count - nrow(dt), "\n",
    "Success rate: ", round(100 * nrow(dt) / initial_count, 2), "%\n"
  )
} else {
  done_content <- paste0(
    done_content,
    "LiftOver: NOT NEEDED (builds match)\n"
  )
}

writeLines(done_content, opt$done)
cat("  Created done marker:", opt$done, "\n")

cat("\n============================================================\n")
cat("LiftOver completed successfully for", opt$dataset, "\n")
cat("============================================================\n\n")
