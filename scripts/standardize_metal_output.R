#!/usr/bin/env Rscript
#
# Standardize METAL meta-analysis output
# - Rename columns to standard format (lowercase for METAL-specific)
# - Split MarkerName back to constituent columns (chrom, pos, snpid)
# - Filter SNPs with p < 5e-3 and create gzipped output
# - Generate Manhattan plot (PDF and PNG) for filtered SNPs
# - Sort by chromosome and position
#

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(ggplot2)
  library(openxlsx)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--input"), type="character", 
              help="Input METAL output file (.tbl)"),
  make_option(c("--output"), type="character",
              help="Output standardized TSV file"),
  make_option(c("--output_filtered"), type="character",
              help="Output filtered gzipped TSV file (p < 5e-3)"),
  make_option(c("--output_excel"), type="character",
              help="Output Excel file for SNPs with p < 5e-5"),
  make_option(c("--plot_pdf"), type="character",
              help="Output Manhattan plot PDF file"),
  make_option(c("--plot_png"), type="character",
              help="Output Manhattan plot PNG file"),
  make_option(c("--combination"), type="character",
              help="Study combination name"),
  make_option(c("--studies"), type="character",
              help="Comma-separated list of study names"),
  make_option(c("--markername_format"), type="character",
              help="MarkerName format pattern (e.g., 'chrom:pos:snpid')"),
  make_option(c("--log"), type="character",
              help="Log file path")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Validate arguments
if (is.null(args$input) || is.null(args$output) || is.null(args$combination)) {
  stop("Required arguments: --input, --output, --combination")
}

# Set up logging
log_file <- args$log
if (!is.null(log_file)) {
  log_con <- file(log_file, open="wt")
  sink(log_con, type="output")
  sink(log_con, type="message")
}

cat("============================================================\n")
cat("Standardizing METAL Output\n")
cat("============================================================\n\n")
cat("Combination:", args$combination, "\n")
cat("Studies:", args$studies, "\n")
cat("Input:", args$input, "\n")
cat("Output:", args$output, "\n")
cat("MarkerName format:", args$markername_format, "\n\n")

# Read METAL output
cat("Reading METAL output...\n")
dt <- fread(args$input, header=TRUE, sep="\t", data.table=TRUE)
cat("  Total variants:", nrow(dt), "\n")
cat("  Input columns:", paste(names(dt), collapse=", "), "\n\n")

# Store total SNP count before filtering
total_snps <- nrow(dt)

# Standardize column names
cat("Standardizing column names...\n")
setnames(dt, old=c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value"),
         new=c("markername", "nea", "ea", "beta", "se", "p"), skip_absent=TRUE)

# Rename METAL-specific columns to lowercase
metal_cols <- c("Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal")
metal_cols_lower <- c("direction", "het_isq", "het_chisq", "het_df", "het_pval")
for (i in seq_along(metal_cols)) {
  if (metal_cols[i] %in% names(dt)) {
    setnames(dt, old=metal_cols[i], new=metal_cols_lower[i])
  }
}

cat("  Standardized columns:", paste(names(dt), collapse=", "), "\n\n")

# Split MarkerName column
cat("Splitting MarkerName column...\n")
if (!is.null(args$markername_format)) {
  format_parts <- strsplit(args$markername_format, ":")[[1]]
  cat("  Format parts:", paste(format_parts, collapse=", "), "\n")
  
  # Split markername by ":"
  markername_split <- tstrsplit(dt$markername, ":", fixed=TRUE, fill=NA, type.convert=FALSE)
  
  # Assign to columns based on format
  for (i in seq_along(format_parts)) {
    col_name <- format_parts[i]
    if (i <= length(markername_split)) {
      dt[, (col_name) := markername_split[[i]]]
    }
  }
  
  # Clean chromosome column - remove "chr" prefix if present
  if ("chrom" %in% names(dt)) {
    dt[, chrom := gsub("^chr", "", chrom, ignore.case=TRUE)]
  }
  
  # Convert position to numeric
  if ("pos" %in% names(dt)) {
    dt[, pos := as.numeric(pos)]
  }
  
  cat("  Created columns:", paste(format_parts, collapse=", "), "\n\n")
}

# Sort by chromosome and position
cat("Sorting by chromosome and position...\n")
if ("chrom" %in% names(dt) && "pos" %in% names(dt)) {
  # Remove scaffold and alternative contig chromosomes (containing "_")
  n_before <- nrow(dt)
  dt_removed <- dt[grepl("_", chrom, fixed=FALSE)]
  n_removed <- nrow(dt_removed)
  
  if (n_removed > 0) {
    cat(paste0("  Removing ", n_removed, " variants from scaffold/alternative contigs:\n"))
    scaffold_summary <- dt_removed[, .N, by=chrom][order(-N)]
    for (i in seq_len(min(10, nrow(scaffold_summary)))) {
      cat(sprintf("    %s: %d variants\n", 
                  scaffold_summary$chrom[i], 
                  scaffold_summary$N[i]))
    }
    if (nrow(scaffold_summary) > 10) {
      cat(sprintf("    ... and %d more scaffold/alt contigs\n", 
                  nrow(scaffold_summary) - 10))
    }
    
    # Remove these variants
    dt <- dt[!grepl("_", chrom, fixed=FALSE)]
    cat(sprintf("  Retained %d variants on standard chromosomes\n\n", nrow(dt)))
  } else {
    cat("  No scaffold/alternative contig chromosomes found\n\n")
  }
  
  # Create numeric chromosome for sorting (1-22, then X, Y, MT)
  dt[, chrom_num := ifelse(
    chrom %in% c("X", "x"), 23,
    ifelse(chrom %in% c("Y", "y"), 24,
    ifelse(chrom %in% c("MT", "mt", "M", "m"), 25, 
    as.numeric(chrom)))
  )]
  
  setorder(dt, chrom_num, pos)
  dt[, chrom_num := NULL]
  
  cat("  Sorted", nrow(dt), "variants\n\n")
}

# Reorder columns: standard first, then METAL-specific
standard_cols <- c("chrom", "pos", "snpid", "ea", "nea", "beta", "se", "p")
metal_specific <- c("direction", "het_isq", "het_chisq", "het_df", "het_pval")
other_cols <- setdiff(names(dt), c(standard_cols, metal_specific, "markername"))

# Build final column order (only include columns that exist)
final_cols <- c(
  intersect(standard_cols, names(dt)),
  intersect(metal_specific, names(dt)),
  other_cols
)
setcolorder(dt, final_cols)

# Write full standardized output
cat("Writing standardized output...\n")
fwrite(dt, args$output, sep="\t", quote=FALSE, na="NA")
cat("  Written:", args$output, "\n")
cat("  Rows:", nrow(dt), "\n")
cat("  Columns:", ncol(dt), "\n\n")

# Filter for significant SNPs (p < 5e-3)
cat("Filtering significant SNPs (p < 5e-3)...\n")
dt_sig <- dt[p < 5e-3]
cat("  Significant variants:", nrow(dt_sig), "\n\n")

# Filter for highly significant SNPs (p < 5e-5) for Excel output
cat("Filtering highly significant SNPs (p < 5e-5)...\n")
dt_excel <- dt[p < 5e-5]
cat("  Highly significant variants:", nrow(dt_excel), "\n\n")

if (nrow(dt_sig) > 0) {
  # Write filtered output (gzipped)
  if (!is.null(args$output_filtered)) {
    cat("Writing filtered output (gzipped)...\n")
    fwrite(dt_sig, args$output_filtered, sep="\t", quote=FALSE, na="NA")
    cat("  Written:", args$output_filtered, "\n\n")
  }
  
  # Write Excel output for highly significant SNPs (p < 5e-5)
  if (!is.null(args$output_excel) && nrow(dt_excel) > 0) {
    cat("Writing Excel output (p < 5e-5)...\n")
    
    # Parse study names
    study_names <- strsplit(args$studies, ",")[[1]]
    study_names <- trimws(study_names)
    cat("  Studies:", paste(study_names, collapse=", "), "\n")
    
    # Prepare Excel data with individual study columns
    excel_dt <- copy(dt_excel)
    
    # Read and merge individual study data
    for (study in study_names) {
      study_file <- paste0("results/04_metal_ready/", study, ".tsv")
      
      if (file.exists(study_file)) {
        cat("  Reading", study, "data...\n")
        study_dt <- fread(study_file, sep="\t", header=TRUE)
        
        # Create merge key
        if ("markername" %in% names(study_dt)) {
          merge_key <- "markername"
        } else if (all(c("chrom", "pos", "snpid") %in% names(study_dt))) {
          study_dt[, markername_tmp := paste(chrom, pos, snpid, sep=":")]
          merge_key <- "markername_tmp"
        } else {
          cat("    WARNING: Cannot merge", study, "- missing key columns\n")
          next
        }
        
        # Select and rename columns for this study
        study_cols <- c("markername", "beta", "se", "p", "eaf")
        if ("or" %in% names(study_dt)) {
          study_cols <- c(study_cols, "or")
        }
        
        # Keep only existing columns
        study_cols <- intersect(study_cols, names(study_dt))
        if (merge_key != "markername") {
          study_cols <- c(merge_key, setdiff(study_cols, "markername"))
        }
        
        study_subset <- study_dt[, ..study_cols]
        
        # Calculate OR from beta if not present
        if (!"or" %in% names(study_subset) && "beta" %in% names(study_subset)) {
          study_subset[, or := exp(beta)]
        }
        
        # Rename columns with study prefix
        old_names <- setdiff(names(study_subset), c("markername", "markername_tmp"))
        new_names <- paste0(study, "_", old_names)
        setnames(study_subset, old_names, new_names)
        
        # Merge with excel_dt
        by_col <- ifelse(merge_key == "markername", "markername", "markername")
        if (merge_key == "markername_tmp") {
          excel_dt[, markername_tmp := markername]
          study_subset[, markername := markername_tmp]
          study_subset[, markername_tmp := NULL]
        }
        
        excel_dt <- merge(excel_dt, study_subset, by="markername", all.x=TRUE)
        
      } else {
        cat("    WARNING: Study file not found:", study_file, "\n")
      }
    }
    
    # Remove temporary columns
    if ("markername_tmp" %in% names(excel_dt)) {
      excel_dt[, markername_tmp := NULL]
    }
    
    # Reorder columns: meta columns first, then individual studies
    meta_cols <- intersect(c("chrom", "pos", "snpid", "ea", "nea", "beta", "se", "p", "eaf",
                              "direction", "het_isq", "het_chisq", "het_df", "het_pval"),
                           names(excel_dt))
    study_cols <- grep("^(MEX|CLM|LAMR|CHI|EUR|MVP)_", names(excel_dt), value=TRUE)
    other_cols <- setdiff(names(excel_dt), c(meta_cols, study_cols))
    
    final_order <- c(meta_cols, study_cols, other_cols)
    setcolorder(excel_dt, final_order)
    
    # Create formatted Excel workbook
    wb <- createWorkbook()
    addWorksheet(wb, "Significant_SNPs")
    
    # Write data
    writeData(wb, "Significant_SNPs", excel_dt, startRow=1, startCol=1)
    
    # Create number formats
    num_2dec <- createStyle(numFmt="0.00")
    sci_format <- createStyle(numFmt="0.00E+00")
    comma_format <- createStyle(numFmt="#,##0")
    
    # Create style for genome-wide significant p-values (light red background)
    sig_pval_style <- createStyle(numFmt="0.00E+00", fgFill="#FFB3B3")
    
    # Apply formatting
    # Position column: comma format
    pos_col <- which(names(excel_dt) == "pos")
    if (length(pos_col) > 0) {
      addStyle(wb, "Significant_SNPs", comma_format, 
               rows=2:(nrow(excel_dt)+1), cols=pos_col, gridExpand=TRUE)
    }
    
    # P-value columns: scientific notation
    p_cols <- which(grepl("^p$|_p$", names(excel_dt)))
    for (col in p_cols) {
      addStyle(wb, "Significant_SNPs", sci_format,
               rows=2:(nrow(excel_dt)+1), cols=col, gridExpand=TRUE)
      
      # Highlight genome-wide significant p-values (p < 5e-8) in light red
      sig_rows <- which(excel_dt[[col]] < 5e-8)
      if (length(sig_rows) > 0) {
        addStyle(wb, "Significant_SNPs", sig_pval_style,
                 rows=sig_rows+1, cols=col, gridExpand=TRUE, stack=FALSE)
      }
    }
    
    # Beta, SE, EAF, OR, het columns: 2 decimals
    num_cols <- which(grepl("^(beta|se|eaf|het_)|_(beta|se|eaf|or)$", names(excel_dt)))
    for (col in num_cols) {
      addStyle(wb, "Significant_SNPs", num_2dec,
               rows=2:(nrow(excel_dt)+1), cols=col, gridExpand=TRUE)
    }
    
    # Set column widths
    setColWidths(wb, "Significant_SNPs", cols=1:ncol(excel_dt), widths="auto")
    
    # Freeze first row
    freezePane(wb, "Significant_SNPs", firstRow=TRUE)
    
    # Save workbook
    saveWorkbook(wb, args$output_excel, overwrite=TRUE)
    cat("  Written:", args$output_excel, "\n")
    cat("  Rows:", nrow(excel_dt), "\n")
    cat("  Columns:", ncol(excel_dt), "\n\n")
    
  } else if (!is.null(args$output_excel) && nrow(dt_excel) == 0) {
    cat("  No SNPs with p < 5e-5, creating empty Excel file\n")
    # Create empty workbook with headers
    wb <- createWorkbook()
    addWorksheet(wb, "Significant_SNPs")
    writeData(wb, "Significant_SNPs", dt[0,])
    saveWorkbook(wb, args$output_excel, overwrite=TRUE)
    cat("  Written:", args$output_excel, "\n\n")
  }
  
  # Create Manhattan plot
  if (!is.null(args$plot_pdf) || !is.null(args$plot_png)) {
    cat("Creating Manhattan plot...\n")
    
    # Prepare data for plotting
    plot_dt <- copy(dt_sig)
    
    # Ensure numeric chromosome for plotting
    plot_dt[, chrom_plot := ifelse(
      chrom %in% c("X", "x"), 23,
      ifelse(chrom %in% c("Y", "y"), 24,
      ifelse(chrom %in% c("MT", "mt", "M", "m"), 25, 
      as.numeric(chrom)))
    )]
    
    # Calculate -log10(p)
    plot_dt[, log10p := -log10(p)]
    
    # Create cumulative position for x-axis
    plot_dt[, chrom_factor := factor(chrom_plot)]
    setorder(plot_dt, chrom_plot, pos)
    
    # Calculate cumulative positions
    chrom_lengths <- plot_dt[, .(max_pos = max(pos, na.rm=TRUE)), by=chrom_plot]
    setorder(chrom_lengths, chrom_plot)
    chrom_lengths[, cumul_start := c(0, cumsum(max_pos)[-.N])]
    plot_dt <- merge(plot_dt, chrom_lengths[, .(chrom_plot, cumul_start)], by="chrom_plot")
    plot_dt[, pos_cumul := pos + cumul_start]
    
    # Calculate center positions for chromosome labels
    axis_dt <- plot_dt[, .(center = mean(pos_cumul)), by=chrom_plot]
    setorder(axis_dt, chrom_plot)
    
    # Create plot title
    plot_title <- sprintf("%s Meta-Analysis\nTotal SNPs: %s | Studies: %s",
                         args$combination,
                         format(total_snps, big.mark=","),
                         args$studies)
    
    # Create Manhattan plot
    p <- ggplot(plot_dt, aes(x=pos_cumul, y=log10p)) +
      geom_point(aes(color=as.factor(chrom_plot)), alpha=0.7, size=1.2) +
      scale_color_manual(values=rep(c("#276FBF", "#183059"), 13)) +
      scale_x_continuous(
        label = axis_dt$chrom_plot,
        breaks = axis_dt$center
      ) +
      geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed", size=0.5) +
      geom_hline(yintercept=-log10(5e-3), color="blue", linetype="dashed", size=0.5) +
      labs(
        title = plot_title,
        x = "Chromosome",
        y = expression(-log[10](P))
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text.x = element_text(angle=0, hjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    # Save as PDF
    if (!is.null(args$plot_pdf)) {
      cat("  Saving PDF:", args$plot_pdf, "\n")
      ggsave(args$plot_pdf, plot=p, width=12, height=6, units="in", dpi=300)
    }
    
    # Save as PNG
    if (!is.null(args$plot_png)) {
      cat("  Saving PNG:", args$plot_png, "\n")
      ggsave(args$plot_png, plot=p, width=12, height=6, units="in", dpi=300)
    }
    
    cat("  Manhattan plot created\n\n")
  }
} else {
  cat("WARNING: No significant SNPs found (p < 5e-3)\n")
  cat("  Skipping filtered output and Manhattan plot\n\n")
  
  # Create empty files to satisfy Snakemake
  if (!is.null(args$output_filtered)) {
    file.create(args$output_filtered)
  }
  if (!is.null(args$plot_pdf)) {
    file.create(args$plot_pdf)
  }
  if (!is.null(args$plot_png)) {
    file.create(args$plot_png)
  }
}

cat("============================================================\n")
cat("METAL output standardization completed\n")
cat("============================================================\n")

# Close log
if (!is.null(log_file)) {
  sink(type="message")
  sink(type="output")
  close(log_con)
}
