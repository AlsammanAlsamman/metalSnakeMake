#!/usr/bin/env Rscript
# Create Excel file from annotated TSV

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

# Check if openxlsx is available and working
use_openxlsx <- FALSE
tryCatch({
  library(openxlsx)
  use_openxlsx <- TRUE
  cat("openxlsx version:", as.character(packageVersion("openxlsx")), "\n")
}, error = function(e) {
  cat("WARNING: openxlsx not available, will use CSV format\n")
})

# Parse arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Input annotated TSV file"),
  make_option(c("--output"), type="character", help="Output Excel file"),
  make_option(c("--pval_threshold"), type="numeric", help="P-value threshold for filtering"),
  make_option(c("--combination"), type="character", help="Combination name")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("=== Creating Excel Output ===\n")
cat("Input:", opt$input, "\n")
cat("Output:", opt$output, "\n")
cat("P-value threshold:", opt$pval_threshold, "\n\n")

# Read annotated TSV
cat("Reading annotated data...\n")
dt <- fread(opt$input, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("Total variants:", nrow(dt), "\n")

# Filter by p-value
cat("Filtering for p <", opt$pval_threshold, "...\n")
excel_dt <- dt[p < opt$pval_threshold]
cat("Variants for Excel:", nrow(excel_dt), "\n\n")

if (nrow(excel_dt) == 0) {
  cat("WARNING: No variants pass p-value threshold. Creating empty Excel.\n")
}

if (!use_openxlsx) {
  # Fallback to CSV
  csv_file <- gsub("\\.xlsx$", ".csv", opt$output)
  cat("Creating CSV file instead:", csv_file, "\n")
  fwrite(excel_dt, csv_file, sep=",", quote=TRUE)
  cat("CSV file created successfully\n")
  cat("\n=== Excel Creation Complete (CSV format) ===\n")
  quit(save="no", status=0)
}

# Create workbook
cat("Creating workbook...\n")
wb <- createWorkbook()
addWorksheet(wb, "ANNOVAR_Annotated")

# Write data
cat("Writing data...\n")
writeData(wb, "ANNOVAR_Annotated", excel_dt)

# Verify data was written to workbook
cat("Verifying data in workbook...\n")
sheet_data <- readWorkbook(wb, sheet = "ANNOVAR_Annotated")
cat("Rows in workbook:", nrow(sheet_data), "\n")
cat("Cols in workbook:", ncol(sheet_data), "\n")

if (nrow(sheet_data) == 0) {
  stop("ERROR: No data was written to workbook!")
}

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
if (nrow(excel_dt) > 0) {
  annotated_rows <- which(!is.na(excel_dt$Func.refGene)) + 1  # +1 for header
  if (length(annotated_rows) > 0) {
    annotated_style <- createStyle(fgFill = "#E6F3FF")
    addStyle(wb, sheet = "ANNOVAR_Annotated", annotated_style, rows = annotated_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE)
    cat("Highlighted", length(annotated_rows), "rows with ANNOVAR annotation\n")
  }
  
  # Highlight genome-wide significant SNPs (p < 5e-8)
  sig_rows <- which(excel_dt$p < 5e-8) + 1
  if (length(sig_rows) > 0) {
    sig_style <- createStyle(fgFill = "#FFE6E6", textDecoration = "bold")
    addStyle(wb, sheet = "ANNOVAR_Annotated", sig_style, rows = sig_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE, stack = TRUE)
    cat("Highlighted", length(sig_rows), "rows with p < 5e-8\n")
  }
}

# Auto-size columns
setColWidths(wb, sheet = "ANNOVAR_Annotated", cols = 1:ncol(excel_dt), widths = "auto")

# Save
cat("Saving Excel file...\n")

# Create a dedicated temp directory in the project (avoid using R default tempdir)
project_temp <- "temp_excel"
dir.create(project_temp, showWarnings = FALSE, recursive = TRUE)
temp_file <- file.path(project_temp, paste0(opt$combination, "_temp.xlsx"))

cat("Writing to temp file:", temp_file, "\n")
saveWorkbook(wb, temp_file, overwrite = TRUE)

# Check temp file
if (file.exists(temp_file)) {
  temp_size <- file.info(temp_file)$size
  cat("Temp file size:", temp_size, "bytes\n")
  
  if (temp_size > 1000) {
    # Copy to final location
    cat("Copying to final location:", opt$output, "\n")
    file.copy(temp_file, opt$output, overwrite = TRUE)
    file.remove(temp_file)
  } else {
    # Debug: try to read the file back
    cat("ERROR: Temp file is only", temp_size, "bytes\n")
    cat("Attempting to diagnose issue...\n")
    
    # Try saving without temp file directly
    cat("Trying direct save...\n")
    saveWorkbook(wb, opt$output, overwrite = TRUE)
    
    final_size <- file.info(opt$output)$size
    cat("Direct save result:", final_size, "bytes\n")
    
    if (final_size < 1000) {
      stop(paste("ERROR: Both temp and direct saves failed. Sizes:", temp_size, "and", final_size, "bytes"))
    }
  }
} else {
  stop("ERROR: Temp file was not created!")
}

# Verify final file
if (file.exists(opt$output)) {
  size <- file.info(opt$output)$size
  cat("SUCCESS! Excel file created:\n")
  cat("  Path:", opt$output, "\n")
  cat("  Size:", size, "bytes\n")
  
  if (size < 1000 && nrow(excel_dt) > 0) {
    stop("ERROR: Excel file is suspiciously small!")
  }
} else {
  stop("ERROR: Excel file was not created!")
}

cat("\n=== Excel Creation Complete ===\n")
