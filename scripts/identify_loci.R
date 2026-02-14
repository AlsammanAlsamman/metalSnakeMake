#!/usr/bin/env Rscript
#
# identify_loci.R
#
# Identifies independent genomic risk loci from GWAS meta-analysis using LD-based clumping.
# 
# Key features:
# 1. Subsets reference panel to SNPs in meta-analysis (cached for efficiency)
# 2. Performs LD clumping to identify independent significant variants  
# 3. Merges correlated signals into genomic loci based on LD and distance
#
# Requirements: R packages data.table, optparse, igraph
#               PLINK accessible as 'plink'

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(igraph)
})

# Parse command line arguments
option_list <- list(
  make_option("--meta", type="character", help="Meta-analysis TSV file (must have: markername, Chr, Pos, P-value, Allele1, Allele2)"),
  make_option("--ref-bfile", type="character", help="PLINK bfile prefix for reference panel"),
  make_option("--ref-name", type="character", help="Reference panel name for cached subset"),
  make_option("--combination", type="character", help="Study combination name"),
  make_option("--temp-dir", type="character", help="Temp directory for reference subset caching"),
  make_option("--outdir", type="character", help="Output directory"),
  make_option("--leadP", type="double", default=5e-8, help="P-value threshold for independent significant SNPs"),
  make_option("--gwasP", type="double", default=1e-5, help="P-value threshold for SNPs in LD window"),
  make_option("--r2", type="double", default=0.6, help="LD r2 threshold for building candidate window"),
  make_option("--r2-2", type="double", default=0.1, help="LD r2 threshold to group independent SNPs into same lead"),
  make_option("--mergeDist", type="integer", default=250, help="Distance (kb) to merge nearby loci"),
  make_option("--maf", type="double", default=0.01, help="Minimum MAF in reference panel"),
  make_option("--refSNPs", type="integer", default=1, help="Include reference-only SNPs (1=yes, 0=no)"),
  make_option("--windowKb", type="integer", default=500, help="Max window (kb) for PLINK LD queries"),
  make_option("--threads", type="integer", default=8, help="Number of threads for PLINK")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate required arguments
if (is.null(opt$meta)) stop("--meta is required")
if (is.null(opt$`ref-bfile`)) stop("--ref-bfile is required")
if (is.null(opt$`ref-name`)) stop("--ref-name is required")
if (is.null(opt$combination)) stop("--combination is required")
if (is.null(opt$`temp-dir`)) stop("--temp-dir is required")
if (is.null(opt$outdir)) stop("--outdir is required")

message("========================================")
message("Genomic Loci Identification Pipeline")
message("========================================")
message(sprintf("Meta-analysis file: %s", opt$meta))
message(sprintf("Reference panel: %s", opt$`ref-name`))
message(sprintf("Combination: %s", opt$combination))
message(sprintf("Lead P threshold: %g", opt$leadP))
message(sprintf("GWAS P threshold: %g", opt$gwasP))
message(sprintf("LD r2 threshold: %g", opt$r2))
message(sprintf("Merge distance: %d kb", opt$mergeDist))

# Create output directories
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$`temp-dir`, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Step 1: Load and prepare meta-analysis data
# =============================================================================
message("\n[Step 1/5] Loading meta-analysis data...")
meta <- fread(opt$meta)
message(sprintf("  Loaded %d variants from meta-analysis", nrow(meta)))

# Normalize column names
colnames(meta) <- trimws(colnames(meta))
col_map <- list(
  markername = c("markername", "snp", "rsid", "id", "variant_id", "MarkerName"),
  chr = c("Chr", "chr", "chromosome", "CHR", "Chromosome"),
  pos = c("Pos", "pos", "position", "bp", "BP", "Position"),
  pval = c("P-value", "p-value", "pval", "p", "P", "Pvalue"),
  a1 = c("Allele1", "allele1", "a1", "effect_allele", "EA"),
  a2 = c("Allele2", "allele2", "a2", "other_allele", "OA")
)

for (target in names(col_map)) {
  candidates <- col_map[[target]]
  hit <- intersect(candidates, colnames(meta))
  if (length(hit) == 0) stop(sprintf("Could not find column for %s in meta file", target))
  if (hit[1] != target) setnames(meta, hit[1], target)
}

# Clean chromosome and position
meta[, chr := gsub("^chr", "", chr, ignore.case = TRUE)]
meta[, chr := gsub("^23$", "X", chr)]
meta[, pos := as.integer(pos)]
meta[, pval := as.numeric(pval)]

# Remove missing values
meta <- meta[!is.na(chr) & !is.na(pos) & !is.na(pval)]
message(sprintf("  Retained %d variants after cleaning", nrow(meta)))

# Sort by chromosome and position
setorder(meta, chr, pos)

# =============================================================================
# Step 2: Subset reference panel to meta SNPs (cached)
# =============================================================================
message("\n[Step 2/5] Subsetting reference panel to meta SNPs...")
subset_prefix <- file.path(opt$`temp-dir`, sprintf("%s_%s_subset", opt$combination, opt$`ref-name`))
subset_bim <- paste0(subset_prefix, ".bim")

if (file.exists(subset_bim)) {
  message(sprintf("  Found cached subset: %s", subset_prefix))
} else {
  message("  Creating reference panel subset (this may take a few minutes)...")
  
  # Write SNP list for extraction
  snp_list_file <- tempfile(pattern = "snp_extract_", fileext = ".txt")
  fwrite(data.table(SNP = unique(meta$markername)), snp_list_file, col.names = FALSE)
  
  # Run PLINK to extract SNPs
  plink_cmd <- sprintf(
    "plink --bfile %s --extract %s --make-bed --out %s --threads %d",
    opt$`ref-bfile`, snp_list_file, subset_prefix, opt$threads
  )
  
  message(sprintf("  Running: %s", plink_cmd))
  plink_result <- system(plink_cmd, intern = FALSE)
  
  if (plink_result != 0) {
    stop("PLINK subset command failed. Check that reference panel exists and is valid.")
  }
  
  unlink(snp_list_file)
  
  if (!file.exists(subset_bim)) {
    stop("Reference panel subset was not created. Check PLINK output above.")
  }
  
  message(sprintf("  Created subset reference panel: %s", subset_prefix))
}

# Load subset BIM to get available SNPs
bim <- fread(subset_bim, header = FALSE)
setnames(bim, c("chr", "rsid", "cm", "bp", "a1", "a2"))
message(sprintf("  Reference subset contains %d SNPs", nrow(bim)))

# Filter meta to SNPs in reference
meta <- meta[markername %chin% bim$rsid]
message(sprintf("  %d meta SNPs matched to reference panel", nrow(meta)))

if (nrow(meta) == 0) stop("No SNPs from meta-analysis found in reference panel")

# =============================================================================
# Step 3: Calculate reference allele frequencies
# =============================================================================
message("\n[Step 3/5] Calculating reference allele frequencies...")
freq_prefix <- file.path(opt$`temp-dir`, sprintf("%s_%s_freq", opt$combination, opt$`ref-name`))
freq_file <- paste0(freq_prefix, ".frq")

if (!file.exists(freq_file)) {
  freq_cmd <- sprintf(
    "plink --bfile %s --freq --out %s --threads %d",
    subset_prefix, freq_prefix, opt$threads
  )
  system(freq_cmd, intern = FALSE)
}

ref_freq <- fread(freq_file)
setnames(ref_freq, toupper(names(ref_freq)))
message(sprintf("  Loaded frequencies for %d SNPs", nrow(ref_freq)))

# =============================================================================
# Step 4: Identify independent significant SNPs
# =============================================================================
message("\n[Step 4/5] Identifying independent significant SNPs...")

# Get SNPs below lead P threshold
sig_snps <- meta[pval < opt$leadP]
setorder(sig_snps, pval, pos)
message(sprintf("  Found %d SNPs with P < %g", nrow(sig_snps), opt$leadP))

if (nrow(sig_snps) == 0) {
  message("  No significant SNPs found. Creating empty output.")
  loci_empty <- data.table(
    GenomicLocus = integer(),
    uniqID = character(),
    rsID = character(),
    chr = character(),
    pos = integer(),
    p = numeric(),
    start = integer(),
    end = integer(),
    nSNPs = integer(),
    nGWASSNPs = integer(),
    nIndSigSNPs = integer(),
    IndSigSNPs = character(),
    leadSNPs = character()
  )
  fwrite(loci_empty, file.path(opt$outdir, "GenomicRiskLoci.txt"), sep = "\t")
  quit(save = "no", status = 0)
}

# Apply MAF filter using reference frequencies
sig_snps_maf <- merge(sig_snps, ref_freq[, .(SNP, MAF)], 
                      by.x = "markername", by.y = "SNP", all.x = TRUE)
sig_snps_maf <- sig_snps_maf[is.na(MAF) | MAF >= opt$maf]
message(sprintf("  %d SNPs pass MAF >= %g filter", nrow(sig_snps_maf), opt$maf))

if (nrow(sig_snps_maf) == 0) {
  message("  No SNPs passed MAF filter. Creating empty output.")
  loci_empty <- data.table(
    GenomicLocus = integer(),
    uniqID = character(),
    rsID = character(),
    chr = character(),
    pos = integer(),
    p = numeric(),
    start = integer(),
    end = integer(),
    nSNPs = integer(),
    nGWASSNPs = integer(),
    nIndSigSNPs = integer(),
    IndSigSNPs = character(),
    leadSNPs = character()
  )
  fwrite(loci_empty, file.path(opt$outdir, "GenomicRiskLoci.txt"), sep = "\t")
  quit(save = "no", status = 0)
}

# Perform LD-based clumping per chromosome
message("  Performing LD-based clumping...")

clump_file <- file.path(opt$`temp-dir`, sprintf("%s_clump_input.txt", opt$combination))
clump_snps <- sig_snps_maf[, .(SNP = markername, P = pval)]
fwrite(clump_snps, clump_file, sep = "\t")

clump_prefix <- file.path(opt$`temp-dir`, sprintf("%s_clump", opt$combination))
clump_cmd <- sprintf(
  "plink --bfile %s --clump %s --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --out %s --threads %d",
  subset_prefix, clump_file, opt$leadP, opt$gwasP, opt$r2, opt$windowKb, clump_prefix, opt$threads
)

system(clump_cmd, intern = FALSE)

clump_result_file <- paste0(clump_prefix, ".clumped")
if (!file.exists(clump_result_file)) {
  warning("PLINK clumping did not produce output. Using all significant SNPs as independent.")
  lead_snps <- sig_snps_maf[, .(markername, chr, pos, pval)]
} else {
  clumped <- fread(clump_result_file)
  lead_snps <- sig_snps_maf[markername %chin% clumped$SNP]
}

message(sprintf("  Identified %d independent lead SNPs", nrow(lead_snps)))

# =============================================================================
# Step 5: Merge leads into genomic loci
# =============================================================================
message("\n[Step 5/5] Merging leads into genomic loci...")

merge_bp <- opt$mergeDist * 1000
loci_list <- list()
locus_id <- 1

# Process each chromosome
chroms <- unique(lead_snps$chr)
for (chrom in chroms) {
  chr_leads <- lead_snps[chr == chrom]
  setorder(chr_leads, pos)
  
  if (nrow(chr_leads) == 0) next
  
  # Merge nearby leads
  current_start <- chr_leads$pos[1]
  current_end <- chr_leads$pos[1]
  current_leads <- chr_leads$markername[1]
  current_top_p <- chr_leads$pval[1]
  current_top_snp <- chr_leads$markername[1]
  
  for (i in seq_len(nrow(chr_leads))) {
    if (i == 1) next
    
    if (chr_leads$pos[i] - current_end <= merge_bp) {
      # Merge into current locus
      current_end <- chr_leads$pos[i]
      current_leads <- paste(current_leads, chr_leads$markername[i], sep = ";")
      if (chr_leads$pval[i] < current_top_p) {
        current_top_p <- chr_leads$pval[i]
        current_top_snp <- chr_leads$markername[i]
      }
    } else {
      # Save current locus
      loci_list[[length(loci_list) + 1]] <- data.table(
        GenomicLocus = locus_id,
        uniqID = sprintf("%s:%d", chrom, current_start),
        rsID = current_top_snp,
        chr = chrom,
        pos = current_start,
        p = current_top_p,
        start = current_start,
        end = current_end,
        nSNPs = length(strsplit(current_leads, ";")[[1]]),
        nGWASSNPs = length(strsplit(current_leads, ";")[[1]]),
        nIndSigSNPs = length(strsplit(current_leads, ";")[[1]]),
        IndSigSNPs = current_leads,
        leadSNPs = current_leads
      )
      locus_id <- locus_id + 1
      
      # Start new locus
      current_start <- chr_leads$pos[i]
      current_end <- chr_leads$pos[i]
      current_leads <- chr_leads$markername[i]
      current_top_p <- chr_leads$pval[i]
      current_top_snp <- chr_leads$markername[i]
    }
  }
  
  # Save final locus for chromosome
  loci_list[[length(loci_list) + 1]] <- data.table(
    GenomicLocus = locus_id,
    uniqID = sprintf("%s:%d", chrom, current_start),
    rsID = current_top_snp,
    chr = chrom,
    pos = current_start,
    p = current_top_p,
    start = current_start,
    end = current_end,
    nSNPs = length(strsplit(current_leads, ";")[[1]]),
    nGWASSNPs = length(strsplit(current_leads, ";")[[1]]),
    nIndSigSNPs = length(strsplit(current_leads, ";")[[1]]),
    IndSigSNPs = current_leads,
    leadSNPs = current_leads
  )
  locus_id <- locus_id + 1
}

# Combine all loci
loci_dt <- rbindlist(loci_list, fill = TRUE)
message(sprintf("  Created %d genomic loci", nrow(loci_dt)))

# Write output
out_file <- file.path(opt$outdir, "GenomicRiskLoci.txt")
fwrite(loci_dt, out_file, sep = "\t", quote = FALSE, na = "NA")
message(sprintf("\n✓ Wrote GenomicRiskLoci.txt to: %s", out_file))

# Summary statistics
message("\n========================================")
message("Summary Statistics")
message("========================================")
message(sprintf("Total genomic loci: %d", nrow(loci_dt)))
message(sprintf("Total independent SNPs: %d", sum(loci_dt$nIndSigSNPs)))
message(sprintf("Median locus size: %d kb", median(loci_dt$end - loci_dt$start) / 1000))
message(sprintf("Most significant P-value: %g", min(loci_dt$p)))
message("========================================")
message("✓ Loci identification completed successfully")
