#!/usr/bin/env Rscript
# Create Excel file from annotated TSV with CSV fallback

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
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

# Try creating Excel file
excel_success <- FALSE

tryCatch({
  library(openxlsx)
  cat("Attempting Excel creation with openxlsx...\n")
  
  wb <- createWorkbook()
  addWorksheet(wb, "ANNOVAR_Annotated")
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
  
  # Highlight variants with annotations
  if (nrow(excel_dt) > 0 && "Func.refGene" %in% names(excel_dt)) {
    annotated_rows <- which(!is.na(excel_dt$Func.refGene)) + 1
    if (length(annotated_rows) > 0) {
      annotated_style <- createStyle(fgFill = "#E6F3FF")
      addStyle(wb, sheet = "ANNOVAR_Annotated", annotated_style, rows = annotated_rows, cols = 1:ncol(excel_dt), gridExpand = TRUE)
      cat("Highlighted", length(annotated_rows), "rows with ANNOVAR annotation\n")
    }
    
    # Highlight genome-wide significant
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
  saveWorkbook(wb, opt$output, overwrite = TRUE)
  
  # Verify
  if (file.exists(opt$output)) {
    size <- file.info(opt$output)$size
    if (size > 1000) {
      cat("SUCCESS! Excel file created:", opt$output, "\n")
      cat("File size:", size, "bytes\n")
      excel_success <- TRUE
    } else {
      cat("WARNING: Excel file is only", size, "bytes (likely corrupt)\n")
    }
  } else {
    cat("WARNING: Excel file was not created\n")
  }
  
}, error = function(e) {
  cat("ERROR creating Excel file:", conditionMessage(e), "\n")
})

# Fallback to CSV if Excel failed
if (!excel_success) {
  cat("\n=== Falling back to CSV format ===\n")
  csv_file <- gsub("\\.xlsx$", ".csv", opt$output)
  cat("Creating CSV file:", csv_file, "\n")
  fwrite(excel_dt, csv_file, sep=",", quote=TRUE)
  
  if (file.exists(csv_file)) {
    size <- file.info(csv_file)$size
    cat("CSV file created successfully:", size, "bytes\n")
    cat("You can open this in Excel.\n")
  }
}

cat("\n=== Output Creation Complete ===\n")
