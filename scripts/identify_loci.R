#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--meta", type = "character", help = "Meta-analysis TSV file"),
  make_option("--ref-bfile", type = "character", help = "PLINK bfile prefix for reference panel"),
  make_option("--ref-name", type = "character", help = "Reference panel name"),
  make_option("--target-build", type = "character", default = "hg38", help = "Target build folder for reference panel (e.g., hg38)"),
  make_option("--combination", type = "character", help = "Study combination name"),
  make_option("--temp-dir", type = "character", help = "Temp directory"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--leadP", type = "double", default = 5e-8),
  make_option("--gwasP", type = "double", default = 1e-5),
  make_option("--r2", type = "double", default = 0.6),
  make_option("--r2-2", type = "double", default = 0.1),
  make_option("--mergeDist", type = "integer", default = 250),
  make_option("--maf", type = "double", default = 0.01),
  make_option("--panel-match-p-threshold", type = "double", default = 5e-3),
  make_option("--refSNPs", type = "integer", default = 1),
  make_option("--windowKb", type = "integer", default = 500),
  make_option("--threads", type = "integer", default = 8),
  make_option("--plink", type = "character", default = "plink")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("meta", "ref-bfile", "ref-name", "combination", "temp-dir", "outdir")
for (name in required) {
  if (is.null(opt[[name]])) stop(sprintf("--%s is required", name))
}

message("========================================")
message("Genomic Loci Identification Pipeline")
message("========================================")
message(sprintf("Meta-analysis file: %s", opt$meta))
message(sprintf("Reference panel: %s", opt$`ref-name`))
message(sprintf("Target build: %s", opt$`target-build`))
message(sprintf("Combination: %s", opt$combination))
message(sprintf("Panel match P threshold: %g", opt$`panel-match-p-threshold`))
message(sprintf("Lead P threshold: %g", opt$leadP))
message(sprintf("GWAS P threshold: %g", opt$gwasP))
message(sprintf("LD r2 threshold: %g", opt$r2))
message(sprintf("Merge distance: %d kb", opt$mergeDist))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$`temp-dir`, showWarnings = FALSE, recursive = TRUE)

col_map <- list(
  markername = c("markername", "MarkerName", "snpid", "snp", "rsid", "id", "variant_id"),
  chr = c("chrom", "CHROM", "Chr", "chr", "chromosome", "CHR", "Chromosome"),
  pos = c("pos", "Pos", "position", "bp", "BP", "Position"),
  pval = c("p", "P", "P-value", "p-value", "pval", "Pvalue"),
  ea = c("ea", "EA", "Allele1", "allele1", "a1", "effect_allele"),
  nea = c("nea", "NEA", "Allele2", "allele2", "a2", "other_allele", "OA")
)

normalize_chr <- function(x) {
  y <- toupper(as.character(x))
  y <- gsub("^CHR", "", y)
  y <- gsub("^23$", "X", y)
  y <- gsub("^24$", "Y", y)
  y
}

resolve_ref_bfile <- function(prefix, target_build) {
  bim <- paste0(prefix, ".bim")
  if (file.exists(bim)) return(prefix)

  prefix_norm <- gsub("\\\\", "/", prefix)
  dir_part <- dirname(prefix_norm)
  base_part <- basename(prefix_norm)
  candidate <- file.path(dir_part, target_build, base_part)
  candidate_bim <- paste0(candidate, ".bim")
  if (file.exists(candidate_bim)) return(candidate)

  stop(sprintf(
    "Reference BIM not found. Checked: %s and %s",
    bim, candidate_bim
  ))
}

write_empty_loci <- function(outdir) {
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
  fwrite(loci_empty, file.path(outdir, "GenomicRiskLoci.txt"), sep = "\t")
}

message("\n[Step 1/6] Loading and normalizing meta-analysis data...")
meta <- fread(opt$meta)
colnames(meta) <- trimws(colnames(meta))

raw_has_snpid <- "snpid" %in% colnames(meta)
if (raw_has_snpid) {
  meta[, snpid_original := as.character(snpid)]
}

for (target in names(col_map)) {
  hit <- intersect(col_map[[target]], colnames(meta))
  if (length(hit) == 0) stop(sprintf("Could not find required column for %s", target))
  if (hit[1] != target) setnames(meta, hit[1], target)
}

meta[, chr := normalize_chr(chr)]
meta[, pos := as.integer(pos)]
meta[, pval := as.numeric(pval)]
meta[, ea := toupper(as.character(ea))]
meta[, nea := toupper(as.character(nea))]
meta[, markername := as.character(markername)]
if (!"snpid" %in% colnames(meta) && "snpid_original" %in% colnames(meta)) {
  meta[, snpid := as.character(snpid_original)]
}
if ("snpid" %in% colnames(meta)) {
  meta[, snpid := as.character(snpid)]
}
meta <- meta[!is.na(chr) & !is.na(pos) & !is.na(pval) & !is.na(ea) & !is.na(nea)]
meta <- meta[pval < opt$`panel-match-p-threshold`]

if (nrow(meta) == 0) {
  message("No variants pass panel match p-value threshold. Writing empty output.")
  write_empty_loci(opt$outdir)
  quit(save = "no", status = 0)
}

if ("snpid" %in% colnames(meta)) {
  meta <- unique(meta[, .(markername, snpid, chr, pos, pval, ea, nea)])
} else {
  meta <- unique(meta[, .(markername, chr, pos, pval, ea, nea)])
}
meta[, meta_id := .I]
setorder(meta, chr, pos, pval)
message(sprintf("  Meta variants retained for panel matching: %d", nrow(meta)))

message("\n[Step 2/6] Matching to reference BIM by chr+pos+alleles...")
ref_bfile_resolved <- resolve_ref_bfile(opt$`ref-bfile`, opt$`target-build`)
message(sprintf("  Resolved reference bfile prefix: %s", ref_bfile_resolved))
bim_file <- paste0(ref_bfile_resolved, ".bim")

bim <- fread(bim_file, header = FALSE)
if (ncol(bim) < 6) stop("Invalid BIM format: expected at least 6 columns")
setnames(bim, c("chr", "panel_id", "cm", "pos", "panel_a1", "panel_a2"))
bim <- bim[, .(chr = normalize_chr(chr), pos = as.integer(pos), panel_id = as.character(panel_id),
               panel_a1 = toupper(as.character(panel_a1)), panel_a2 = toupper(as.character(panel_a2)))]

cand <- merge(meta, bim, by = c("chr", "pos"), all = FALSE, allow.cartesian = TRUE)
if (nrow(cand) == 0) {
  stop("No chr+pos overlaps between meta variants and reference BIM")
}

cand[, direct_match := (ea == panel_a1 & nea == panel_a2)]
cand[, swapped_match := (ea == panel_a2 & nea == panel_a1)]
cand <- cand[direct_match | swapped_match]

if (nrow(cand) == 0) {
  stop("No allele-compatible matches found (direct or swapped) between meta and reference BIM")
}

cand[, match_type := ifelse(direct_match, "direct", "swapped")]
setorder(cand, meta_id, -direct_match)

ambiguous_by_meta <- cand[, .N, by = meta_id][N > 1]
selected <- cand[, .SD[1], by = meta_id]
selected[, markernamepanel := panel_id]

overlap_count <- nrow(merge(meta, unique(bim[, .(chr, pos)]), by = c("chr", "pos"), all = FALSE))
match_summary <- data.table(
  metric = c(
    "meta_variants_tested",
    "chr_pos_overlap_rows",
    "allele_compatible_rows",
    "direct_match_rows",
    "swapped_match_rows",
    "ambiguous_meta_variants",
    "unique_meta_variants_matched",
    "match_rate_percent"
  ),
  value = c(
    nrow(meta),
    overlap_count,
    nrow(cand),
    cand[match_type == "direct", .N],
    cand[match_type == "swapped", .N],
    nrow(ambiguous_by_meta),
    nrow(selected),
    round(100 * nrow(selected) / nrow(meta), 2)
  )
)

message(sprintf("  Matched unique variants: %d / %d (%.2f%%)",
                nrow(selected), nrow(meta), 100 * nrow(selected) / nrow(meta)))
message(sprintf("  Direct matches: %d | Swapped matches: %d | Ambiguous meta variants: %d",
                cand[match_type == "direct", .N], cand[match_type == "swapped", .N], nrow(ambiguous_by_meta)))

if ("snpid" %in% colnames(selected)) {
  panel_match_table <- selected[, .(
    chr,
    pos,
    snpid,
    markername,
    markernamepanel,
    ea,
    nea,
    panel_a1,
    panel_a2,
    match_type,
    pval
  )]
} else {
  panel_match_table <- selected[, .(
    chr,
    pos,
    markername,
    markernamepanel,
    ea,
    nea,
    panel_a1,
    panel_a2,
    match_type,
    pval
  )]
}
setorder(panel_match_table, pval, chr, pos)

panel_match_file <- file.path(opt$outdir, sprintf("%s_panel_markername.tsv", opt$combination))
summary_file <- file.path(opt$outdir, sprintf("%s_panel_markername_summary.tsv", opt$combination))
fwrite(panel_match_table, panel_match_file, sep = "\t", quote = FALSE)
fwrite(match_summary, summary_file, sep = "\t", quote = FALSE)

message(sprintf("  Wrote panel marker table: %s", panel_match_file))
message(sprintf("  Wrote panel match summary: %s", summary_file))

message("\n[Step 3/6] Building reference subset by markernamepanel...")
subset_prefix <- file.path(opt$`temp-dir`, sprintf("%s_%s_panel_subset", opt$combination, opt$`ref-name`))
subset_bim <- paste0(subset_prefix, ".bim")
extract_list <- file.path(opt$`temp-dir`, sprintf("%s_panel_extract_ids.txt", opt$combination))
fwrite(unique(panel_match_table[, .(SNP = markernamepanel)]), extract_list, col.names = FALSE)

plink_subset_cmd <- sprintf(
  "%s --bfile %s --extract %s --make-bed --allow-extra-chr --out %s --threads %d",
  opt$plink, ref_bfile_resolved, extract_list, subset_prefix, opt$threads
)

message(sprintf("  Running: %s", plink_subset_cmd))
subset_status <- system(plink_subset_cmd, intern = FALSE)
if (subset_status != 0 || !file.exists(subset_bim)) {
  stop("PLINK subset creation failed for markernamepanel IDs")
}

bim_subset <- fread(subset_bim, header = FALSE)
setnames(bim_subset, c("chr", "panel_id", "cm", "pos", "a1", "a2"))
panel_match_table <- panel_match_table[markernamepanel %chin% bim_subset$panel_id]

if (nrow(panel_match_table) == 0) {
  stop("No matched markernamepanel IDs found in subset BIM after extraction")
}

message(sprintf("  Markernamepanel IDs retained after subset: %d", nrow(panel_match_table)))

message("\n[Step 4/6] Calculating MAF and selecting significant SNPs...")
freq_prefix <- file.path(opt$`temp-dir`, sprintf("%s_%s_panel_freq", opt$combination, opt$`ref-name`))
freq_file <- paste0(freq_prefix, ".frq")

if (!file.exists(freq_file)) {
  freq_cmd <- sprintf(
    "%s --bfile %s --freq --allow-extra-chr --out %s --threads %d",
    opt$plink, subset_prefix, freq_prefix, opt$threads
  )
  system(freq_cmd, intern = FALSE)
}

if (!file.exists(freq_file)) stop("PLINK frequency calculation failed")

ref_freq <- fread(freq_file)
setnames(ref_freq, toupper(names(ref_freq)))

meta_panel <- merge(
  panel_match_table,
  ref_freq[, .(SNP, MAF)],
  by.x = "markernamepanel",
  by.y = "SNP",
  all.x = TRUE
)
meta_panel <- meta_panel[is.na(MAF) | MAF >= opt$maf]

sig_snps <- meta_panel[pval < opt$leadP]
setorder(sig_snps, pval, chr, pos)

if (nrow(sig_snps) == 0) {
  message("No SNPs pass leadP after panel matching and MAF filter. Writing empty output.")
  write_empty_loci(opt$outdir)
  quit(save = "no", status = 0)
}

message(sprintf("  Significant SNPs for clumping: %d", nrow(sig_snps)))

message("\n[Step 5/6] PLINK clumping using markernamepanel...")
clump_file <- file.path(opt$`temp-dir`, sprintf("%s_clump_input_panel.txt", opt$combination))
clump_prefix <- file.path(opt$`temp-dir`, sprintf("%s_clump_panel", opt$combination))
fwrite(unique(meta_panel[, .(SNP = markernamepanel, P = pval)]), clump_file, sep = "\t")

clump_cmd <- sprintf(
  "%s --bfile %s --clump %s --clump-p1 %g --clump-p2 %g --clump-r2 %g --clump-kb %d --allow-extra-chr --out %s --threads %d",
  opt$plink, subset_prefix, clump_file, opt$leadP, opt$gwasP, opt$r2, opt$windowKb, clump_prefix, opt$threads
)
system(clump_cmd, intern = FALSE)

clump_result_file <- paste0(clump_prefix, ".clumped")
if (!file.exists(clump_result_file)) {
  warning("PLINK clumping output missing; fallback to all significant SNPs")
  lead_snps <- sig_snps
} else {
  clumped <- fread(clump_result_file)
  lead_snps <- sig_snps[markernamepanel %chin% clumped$SNP]
  if (nrow(lead_snps) == 0) {
    warning("No lead SNPs from clump matched back; fallback to all significant SNPs")
    lead_snps <- sig_snps
  }
}

message(sprintf("  Independent lead SNPs: %d", nrow(lead_snps)))

if ("snpid" %in% colnames(lead_snps)) {
  lead_snps[, lead_id := fifelse(!is.na(snpid) & snpid != "", snpid, markernamepanel)]
} else {
  lead_snps[, lead_id := markernamepanel]
}

message("\n[Step 6/6] Merging lead SNPs into genomic loci...")
merge_bp <- opt$mergeDist * 1000
setorder(lead_snps, chr, pos, pval)

loci_list <- list()
locus_id <- 1

for (chrom in unique(lead_snps$chr)) {
  chr_leads <- lead_snps[chr == chrom]
  if (nrow(chr_leads) == 0) next

  current_start <- chr_leads$pos[1]
  current_end <- chr_leads$pos[1]
  current_top_p <- chr_leads$pval[1]
  current_top_snp <- chr_leads$lead_id[1]
  current_leads <- chr_leads$lead_id[1]

  if (nrow(chr_leads) > 1) {
    for (i in 2:nrow(chr_leads)) {
      if ((chr_leads$pos[i] - current_end) <= merge_bp) {
        current_end <- chr_leads$pos[i]
        current_leads <- paste(current_leads, chr_leads$lead_id[i], sep = ";")
        if (chr_leads$pval[i] < current_top_p) {
          current_top_p <- chr_leads$pval[i]
          current_top_snp <- chr_leads$lead_id[i]
        }
      } else {
        loci_list[[length(loci_list) + 1]] <- data.table(
          GenomicLocus = locus_id,
          uniqID = sprintf("%s:%d", chrom, current_start),
          rsID = current_top_snp,
          chr = chrom,
          pos = current_start,
          p = current_top_p,
          start = current_start,
          end = current_end,
          nSNPs = length(strsplit(current_leads, ";", fixed = TRUE)[[1]]),
          nGWASSNPs = length(strsplit(current_leads, ";", fixed = TRUE)[[1]]),
          nIndSigSNPs = length(strsplit(current_leads, ";", fixed = TRUE)[[1]]),
          IndSigSNPs = current_leads,
          leadSNPs = current_leads
        )
        locus_id <- locus_id + 1

        current_start <- chr_leads$pos[i]
        current_end <- chr_leads$pos[i]
        current_top_p <- chr_leads$pval[i]
        current_top_snp <- chr_leads$lead_id[i]
        current_leads <- chr_leads$lead_id[i]
      }
    }
  }

  loci_list[[length(loci_list) + 1]] <- data.table(
    GenomicLocus = locus_id,
    uniqID = sprintf("%s:%d", chrom, current_start),
    rsID = current_top_snp,
    chr = chrom,
    pos = current_start,
    p = current_top_p,
    start = current_start,
    end = current_end,
    nSNPs = length(strsplit(current_leads, ";", fixed = TRUE)[[1]]),
    nGWASSNPs = length(strsplit(current_leads, ";", fixed = TRUE)[[1]]),
    nIndSigSNPs = length(strsplit(current_leads, ";", fixed = TRUE)[[1]]),
    IndSigSNPs = current_leads,
    leadSNPs = current_leads
  )
  locus_id <- locus_id + 1
}

loci_dt <- if (length(loci_list) > 0) rbindlist(loci_list, fill = TRUE) else data.table(
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

out_file <- file.path(opt$outdir, "GenomicRiskLoci.txt")
fwrite(loci_dt, out_file, sep = "\t", quote = FALSE, na = "NA")

message("========================================")
message("Summary Statistics")
message("========================================")
message(sprintf("Total genomic loci: %d", nrow(loci_dt)))
if (nrow(loci_dt) > 0) {
  message(sprintf("Total independent SNPs: %d", sum(loci_dt$nIndSigSNPs)))
  message(sprintf("Median locus size: %.2f kb", median(loci_dt$end - loci_dt$start) / 1000))
  message(sprintf("Most significant P-value: %g", min(loci_dt$p)))
}
message(sprintf("Panel marker table: %s", panel_match_file))
message(sprintf("Panel marker summary: %s", summary_file))
message(sprintf("Genomic loci output: %s", out_file))
message("✓ Loci identification completed successfully")
