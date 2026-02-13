#!/usr/bin/env Rscript
# Calculate Standard Error (SE) from beta, eaf, n_cases, n_controls
# If SE already exists, copy file unchanged

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--dataset"), type="character", help="Dataset name"),
  make_option(c("--input"), type="character", help="Input TSV file"),
  make_option(c("--output"), type="character", help="Output TSV file"),
  make_option(c("--done"), type="character", help="Done marker file"),
  make_option(c("--log"), type="character", help="Log file"),
  make_option(c("--n_cases"), type="integer", help="Number of cases"),
  make_option(c("--n_controls"), type="integer", help="Number of controls")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to calculate SE from beta and p-value
calc_se_from_beta_p <- function(beta, p) {
  # Convert p-value to Z (two-sided)
  Z <- qnorm(1 - p/2)
  
  # Restore sign
  Z <- sign(beta) * Z
  
  # Standard error
  SE <- abs(beta / Z)
  
  return(SE)
}

# Main processing
cat("\n============================================================\n")
cat("Calculate Missing Columns:", opt$dataset, "\n")
cat("============================================================\n\n")

cat("Input file:", opt$input, "\n")
cat("Sample size:", opt$n_cases, "cases +", opt$n_controls, "controls =", 
    opt$n_cases + opt$n_controls, "total\n\n")

# Read input file
cat("Reading input file...\n")
dt <- fread(opt$input, sep="\t", header=TRUE)
cat("  Loaded", nrow(dt), "variants with", ncol(dt), "columns\n")

# Check if beta column exists, if not try to create from OR
has_beta <- "beta" %in% colnames(dt)
has_or <- "or" %in% colnames(dt)
beta_calculated_from_or <- FALSE

if (!has_beta && has_or) {
  cat("\n  Beta column missing but OR column present\n")
  cat("  Calculating beta = log(OR)...\n")
  
  # Calculate beta from OR
  valid_or <- !is.na(dt$or) & dt$or > 0
  n_valid_or <- sum(valid_or)
  
  if (n_valid_or > 0) {
    dt[valid_or, beta := log(or)]
    cat("    Calculated beta for", n_valid_or, "variants with valid OR > 0\n")
    has_beta <- TRUE
    beta_calculated_from_or <- TRUE
  } else {
    cat("    WARNING: No valid OR values (need OR > 0)\n")
  }
} else if (has_beta) {
  cat("\n  Beta column is present\n")
} else {
  cat("\n  WARNING: Neither beta nor OR column found\n")
}

# Check if SE column exists and is not all NA
has_se <- "se" %in% colnames(dt)
se_missing <- FALSE
se_calculated <- FALSE

if (has_se) {
  # Check if SE column has values
  n_non_na <- sum(!is.na(dt$se))
  if (n_non_na == 0) {
    cat("\n  SE column exists but is completely NA (", nrow(dt), "missing values)\n")
    se_missing <- TRUE
  } else {
    cat("\n  SE column exists with", n_non_na, "non-NA values\n")
    cat("  No calculation needed - copying file unchanged\n")
  }
} else {
  cat("\n  SE column is missing from dataset\n")
  se_missing <- TRUE
}

# Calculate SE if missing
if (se_missing) {
  cat("\nCalculating SE from beta and p-value...\n")
  
  # Check required columns
  required_cols <- c("beta", "p")
  missing_req <- setdiff(required_cols, colnames(dt))
  
  if (length(missing_req) > 0) {
    stop("ERROR: Cannot calculate SE - missing required columns: ", 
         paste(missing_req, collapse=", "))
  }
  
  # Calculate SE for non-NA beta and p values
  valid_idx <- !is.na(dt$beta) & !is.na(dt$p) & dt$p > 0 & dt$p < 1
  n_valid <- sum(valid_idx)
  
  if (n_valid == 0) {
    stop("ERROR: No valid variants for SE calculation (need non-NA beta, 0 < p < 1)")
  }
  
  cat("  Valid variants for calculation:", n_valid, "/", nrow(dt), "\n")
  
  # Calculate SE
  dt[valid_idx, se := calc_se_from_beta_p(beta, p)]
  
  n_calculated <- sum(!is.na(dt$se))
  cat("  Calculated SE for", n_calculated, "variants\n")
  
  se_calculated <- TRUE
  
  # Generate density plot for calculated SE
  cat("\nGenerating density plot for calculated SE...\n")
  plot_file <- sub("\\.tsv$", "_se_density.pdf", opt$output)
  
  se_values <- dt[!is.na(se), se]
  
  pdf(plot_file, width=10, height=7)
  p <- ggplot(data.frame(se=se_values), aes(x=se)) +
    geom_density(fill="steelblue", alpha=0.6) +
    geom_vline(xintercept=median(se_values), linetype="dashed", color="red") +
    labs(
      title=paste("Density Plot of Calculated SE -", opt$dataset),
      subtitle=paste("n =", length(se_values), "variants | Median SE =", 
                     round(median(se_values), 4)),
      x="Standard Error (SE)",
      y="Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size=14, face="bold"),
      plot.subtitle = element_text(size=10)
    )
  print(p)
  dev.off()
  
  cat("  Saved density plot to:", plot_file, "\n")
  cat("  SE summary statistics:\n")
  cat("    Min    :", round(min(se_values), 4), "\n")
  cat("    Q1     :", round(quantile(se_values, 0.25), 4), "\n")
  cat("    Median :", round(median(se_values), 4), "\n")
  cat("    Mean   :", round(mean(se_values), 4), "\n")
  cat("    Q3     :", round(quantile(se_values, 0.75), 4), "\n")
  cat("    Max    :", round(max(se_values), 4), "\n")
}

# Write output
cat("\nWriting output file...\n")
dir.create(dirname(opt$output), recursive=TRUE, showWarnings=FALSE)
fwrite(dt, opt$output, sep="\t", quote=FALSE)
cat("  Wrote", nrow(dt), "variants to:", opt$output, "\n")

# Write .done marker
done_content <- paste0(
  "Calculation completed for ", opt$dataset, "\n",
  "Input: ", opt$input, "\n",
  "Output: ", opt$output, "\n",
  "Variants: ", nrow(dt), "\n",
  "Columns: ", paste(colnames(dt), collapse=", "), "\n"
)

if (se_calculated) {
  done_content <- paste0(done_content, "SE_calculated: TRUE\n")
  done_content <- paste0(done_content, "SE_median: ", round(median(dt$se, na.rm=TRUE), 4), "\n")
  if (beta_calculated_from_or) {
    done_content <- paste0(done_content, "Beta_calculated_from_OR: TRUE\n")
  }
} else {
  done_content <- paste0(done_content, "SE_calculated: FALSE (already present)\n")
}

writeLines(done_content, opt$done)
cat("  Created done marker:", opt$done, "\n")

cat("\n============================================================\n")
cat("Calculation completed successfully for", opt$dataset, "\n")
cat("============================================================\n\n")
