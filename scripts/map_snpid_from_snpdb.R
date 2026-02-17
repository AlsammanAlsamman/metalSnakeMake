#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(optparse)
})

option_list <- list(
  make_option("--dataset", type = "character"),
  make_option("--build", type = "character"),
  make_option("--input", type = "character"),
  make_option("--output", type = "character"),
  make_option("--done", type = "character"),
  make_option("--summary", type = "character"),
  make_option("--log", type = "character"),
  make_option("--snpdb-root", type = "character", default = "resources/SNPdb"),
  make_option("--temp-dir", type = "character"),
  make_option("--threads", type = "integer", default = 2),
  make_option("--unmapped-pos-limit", type = "integer", default = 50)
)

opt <- parse_args(OptionParser(option_list = option_list))

get_opt <- function(opts, keys) {
  for (k in keys) {
    val <- opts[[k]]
    if (!is.null(val) && (!is.character(val) || nzchar(val))) {
      return(val)
    }
  }
  NULL
}

# Normalize option aliases for hyphen/dot/underscore naming differences in optparse
opt$snpdb_root <- get_opt(opt, c("snpdb_root", "snpdb-root", "snpdb.root"))
opt$temp_dir <- get_opt(opt, c("temp_dir", "temp-dir", "temp.dir"))

unmapped_file <- file.path(dirname(opt$summary), sprintf("%s_unmapped.tsv", opt$dataset))
unmapped_pos_lookup_file <- file.path(dirname(opt$summary), sprintf("%s_unmapped_pos_lookup.tsv", opt$dataset))

required <- c("dataset", "build", "input", "output", "done", "summary", "log", "snpdb_root", "temp_dir")
missing_required <- required[sapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]))]
if (length(missing_required) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_required, collapse = ", ")))
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$log), recursive = TRUE, showWarnings = FALSE)
dir.create(opt$temp_dir, recursive = TRUE, showWarnings = FALSE)

normalize_chr <- function(x) {
  x <- toupper(as.character(x))
  x <- sub("^CHR", "", x)
  x[x == "M"] <- "MT"
  x
}

complement_allele <- function(x) {
  x <- toupper(as.character(x))
  chartr("ACGT", "TGCA", x)
}

build_to_keys <- function(build) {
  b <- tolower(build)
  if (b %in% c("hg19", "grch37", "37")) {
    return(list(by_chr = "GRCh37", split = "hg19", aliases = c("hg19", "GRCh37")))
  }
  if (b %in% c("hg38", "grch38", "38")) {
    return(list(by_chr = "GRCh38", split = "hg38", aliases = c("hg38", "GRCh38")))
  }
  stop(sprintf("Unsupported build '%s'. Expected hg19/GRCh37 or hg38/GRCh38.", build))
}

find_chr_vcf <- function(snpdb_root, build_info, chrom_norm) {
  chrom_candidates <- unique(c(
    chrom_norm,
    paste0("chr", chrom_norm),
    if (chrom_norm == "MT") c("M", "chrM") else character(0)
  ))

  path_candidates <- character(0)

  for (chr_name in chrom_candidates) {
    path_candidates <- c(
      path_candidates,
      file.path(snpdb_root, "by_chr", build_info$by_chr, paste0(chr_name, ".vcf.gz")),
      file.path(snpdb_root, build_info$split, paste0(chr_name, ".vcf.gz")),
      file.path(snpdb_root, build_info$aliases[1], paste0(chr_name, ".vcf.gz")),
      file.path(snpdb_root, build_info$aliases[2], paste0(chr_name, ".vcf.gz"))
    )
  }

  path_candidates <- unique(path_candidates)
  hits <- path_candidates[file.exists(path_candidates)]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

query_db_for_positions <- function(vcf_file, chrom_norm, pos_dt, temp_regions_file) {
  contig_candidates <- unique(c(
    chrom_norm,
    paste0("chr", chrom_norm),
    if (chrom_norm == "MT") c("M", "chrM") else character(0)
  ))

  for (contig in contig_candidates) {
    reg <- copy(pos_dt)
    reg[, chrom := contig]
    setcolorder(reg, c("chrom", "pos", "end"))
    fwrite(reg, temp_regions_file, sep = "\t", col.names = FALSE)

    cmd <- sprintf(
      "bcftools query -R %s -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\n' %s",
      shQuote(temp_regions_file),
      shQuote(vcf_file)
    )

    db <- tryCatch(
      fread(cmd = cmd, sep = "\t", header = FALSE, showProgress = FALSE),
      error = function(e) data.table()
    )

    if (nrow(db) > 0) {
      setnames(db, c("db_chrom", "pos", "rsid", "ref", "alt"))
      return(db)
    }
  }

  data.table(db_chrom = character(), pos = integer(), rsid = character(), ref = character(), alt = character())
}

map_one_chrom <- function(chr_norm, dt_chr, snpdb_root, build_info, temp_dir) {
  chr_dir <- file.path(temp_dir, paste0("chr_", chr_norm))
  dir.create(chr_dir, recursive = TRUE, showWarnings = FALSE)

  chunk_file <- file.path(chr_dir, "input_chunk.tsv.gz")
  pos_file <- file.path(chr_dir, "positions.tsv")
  db_file <- file.path(chr_dir, "snpdb_subset.tsv.gz")

  fwrite(dt_chr, chunk_file, sep = "\t")

  pos_dt <- unique(dt_chr[, .(pos = as.integer(pos))])
  pos_dt <- pos_dt[!is.na(pos)]
  pos_dt[, end := pos]
  if (nrow(pos_dt) == 0) {
    dt_chr[, mapped_rsid := NA_character_]
    dt_chr[, map_status := "no_pos"]
    return(dt_chr)
  }

  vcf_file <- find_chr_vcf(snpdb_root, build_info, chr_norm)
  if (is.na(vcf_file)) {
    dt_chr[, mapped_rsid := NA_character_]
    dt_chr[, map_status := "missing_snpdb_chr"]
    return(dt_chr)
  }

  db <- query_db_for_positions(vcf_file, chr_norm, pos_dt, pos_file)

  if (nrow(db) == 0) {
    dt_chr[, mapped_rsid := NA_character_]
    dt_chr[, map_status := "no_db_match_pos"]
    return(dt_chr)
  }

  fwrite(db, db_file, sep = "\t")

  db <- db[!is.na(rsid) & rsid != "."]
  if (nrow(db) == 0) {
    dt_chr[, mapped_rsid := NA_character_]
    dt_chr[, map_status := "no_valid_rsid"]
    return(dt_chr)
  }

  db <- db[, .(pos = as.integer(pos), rsid, ref = toupper(ref), alt = toupper(alt))]
  db <- db[!grepl(",", alt)]

  dt_chr[, `:=`(
    ea_u = toupper(as.character(ea)),
    nea_u = toupper(as.character(nea)),
    ea_c = complement_allele(ea),
    nea_c = complement_allele(nea)
  )]

  setkey(dt_chr, pos)
  setkey(db, pos)
  cands <- db[dt_chr, allow.cartesian = TRUE]

  cands[, rank := fifelse(ea_u == alt & nea_u == ref, 1L,
                   fifelse(ea_u == ref & nea_u == alt, 2L,
                   fifelse(ea_c == alt & nea_c == ref, 3L,
                   fifelse(ea_c == ref & nea_c == alt, 4L, 99L))))]

  cands <- cands[rank < 99L]

  if (nrow(cands) == 0) {
    dt_chr[, `:=`(mapped_rsid = NA_character_, map_status = "no_allele_match")]
    dt_chr[, c("ea_u", "nea_u", "ea_c", "nea_c") := NULL]
    return(dt_chr)
  }

  best <- cands[order(row_id, rank)][, .SD[1], by = row_id]
  best <- best[, .(row_id, mapped_rsid = rsid, map_rank = rank)]

  setkey(dt_chr, row_id)
  setkey(best, row_id)
  dt_chr <- best[dt_chr]

  dt_chr[, map_status := fifelse(is.na(mapped_rsid), "unmapped",
                          fifelse(map_rank == 1L, "direct",
                          fifelse(map_rank == 2L, "swapped",
                          fifelse(map_rank == 3L, "strand_direct", "strand_swapped"))))]

  dt_chr[, c("ea_u", "nea_u", "ea_c", "nea_c", "map_rank") := NULL]
  dt_chr
}

cat(sprintf("\n============================================================\n"))
cat(sprintf("SNPdb rsID mapping for dataset: %s\n", opt$dataset))
cat(sprintf("============================================================\n"))
cat(sprintf("Build: %s\n", opt$build))
cat(sprintf("Input: %s\nOutput: %s\n", opt$input, opt$output))
cat(sprintf("SNPdb root: %s\n", opt$snpdb_root))
cat(sprintf("Temp dir: %s\nThreads: %d\n\n", opt$temp_dir, opt$threads))

dt <- fread(opt$input, sep = "\t", showProgress = TRUE)

# If dataset already has snpid column, copy directly as requested.
if ("snpid" %in% names(dt)) {
  total_n <- nrow(dt)
  non_missing_snpid <- dt[!is.na(snpid) & trimws(as.character(snpid)) != "", .N]

  fwrite(dt, opt$output, sep = "\t", na = "NA")

  empty_unmapped <- dt[0]
  empty_unmapped[, `:=`(mapped_rsid = character(), map_status = character())]
  fwrite(empty_unmapped, unmapped_file, sep = "\t", na = "NA")

  empty_lookup <- data.table(
    row_id = integer(),
    chrom = character(),
    pos = integer(),
    ea = character(),
    nea = character(),
    map_status = character(),
    pos_only_rsid = character()
  )
  fwrite(empty_lookup, unmapped_pos_lookup_file, sep = "\t", na = "NA")

  summary_dt <- data.table(
    dataset = opt$dataset,
    build = opt$build,
    total_variants = total_n,
    mapped_rsid_variants = 0L,
    existing_snpid_variants = non_missing_snpid,
    mode = "copied_existing_snpid_column"
  )
  fwrite(summary_dt, opt$summary, sep = "\t")

  with(file(opt$done, "w"), {
    writeLines(c(
      sprintf("SNPdb mapping skipped (existing snpid column) for %s", opt$dataset),
      sprintf("Build: %s", opt$build),
      sprintf("Input: %s", opt$input),
      sprintf("Output: %s", opt$output),
      sprintf("Summary: %s", opt$summary),
      sprintf("Unmapped list: %s", unmapped_file),
      sprintf("Unmapped chr+pos lookup list: %s", unmapped_pos_lookup_file),
      sprintf("Total variants: %d", total_n),
      sprintf("Unmapped variants: %d", 0),
      sprintf("Existing snpid variants (non-missing): %d", non_missing_snpid)
    ))
  })

  cat(sprintf("Dataset %s already contains snpid column; copied without remapping.\n", opt$dataset))
  cat(sprintf("Existing non-missing snpid: %d / %d\n", non_missing_snpid, total_n))
  cat("Done.\n")
  quit(save = "no", status = 0)
}

required_cols <- c("chrom", "pos", "ea", "nea")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop(sprintf("Input missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

if (!"snpid" %in% names(dt)) {
  dt[, snpid := NA_character_]
}

dt[, row_id := .I]
dt[, chrom_norm := normalize_chr(chrom)]
dt[, pos := suppressWarnings(as.integer(pos))]

chroms <- unique(dt[!is.na(chrom_norm), chrom_norm])
chroms <- chroms[!is.na(chroms) & nzchar(chroms)]

if (length(chroms) == 0) {
  stop("No valid chromosomes found in input.")
}

build_info <- build_to_keys(opt$build)
workers <- max(1L, min(as.integer(opt$threads), length(chroms)))
unmapped_pos_limit <- max(0L, as.integer(opt$`unmapped-pos-limit`))

chr_results <- mclapply(
  chroms,
  function(chr_norm) {
    chr_dt <- dt[chrom_norm == chr_norm]
    map_one_chrom(chr_norm, chr_dt, opt$snpdb_root, build_info, opt$temp_dir)
  },
  mc.cores = workers
)

mapped_dt <- rbindlist(chr_results, use.names = TRUE, fill = TRUE)
setorder(mapped_dt, row_id)
out <- mapped_dt

# Keep original alleles as requested; update snpid after all mapping passes.
out[, snpid_orig := snpid]

# Pass 2: rescue only currently-unmapped variants using multiallelic ALT splitting.
pre_rescue_unmapped <- out[is.na(mapped_rsid) | mapped_rsid == ""]
rescued_count <- 0L

if (nrow(pre_rescue_unmapped) > 0) {
  rescue_chunks <- list()
  rescue_chr <- unique(pre_rescue_unmapped[!is.na(chrom_norm) & !is.na(pos), chrom_norm])

  for (chr_norm in rescue_chr) {
    chr_rows <- pre_rescue_unmapped[chrom_norm == chr_norm & !is.na(pos)]
    if (nrow(chr_rows) == 0) next

    vcf_file <- find_chr_vcf(opt$snpdb_root, build_info, chr_norm)
    if (is.na(vcf_file)) next

    pos_dt <- unique(chr_rows[, .(pos = as.integer(pos))])
    pos_dt <- pos_dt[!is.na(pos)]
    if (nrow(pos_dt) == 0) next
    pos_dt[, end := pos]

    reg_file <- file.path(opt$temp_dir, sprintf("multiallelic_rescue_%s.tsv", chr_norm))
    db <- query_db_for_positions(vcf_file, chr_norm, pos_dt, reg_file)
    if (nrow(db) == 0) next

    db <- db[!is.na(rsid) & rsid != ".", .(pos = as.integer(pos), rsid, ref = toupper(ref), alt = toupper(alt))]
    db <- db[!is.na(ref) & ref != "" & !is.na(alt) & alt != ""]
    if (nrow(db) == 0) next

    db_exp <- db[, .(alt_allele = unlist(strsplit(alt, ",", fixed = TRUE))), by = .(pos, rsid, ref)]
    db_exp <- db_exp[!is.na(alt_allele) & alt_allele != "" & alt_allele != "."]
    if (nrow(db_exp) == 0) next

    chr_rows[, `:=`(
      ea_u = toupper(as.character(ea)),
      nea_u = toupper(as.character(nea)),
      ea_c = complement_allele(ea),
      nea_c = complement_allele(nea)
    )]

    cand <- merge(
      chr_rows[, .(row_id, pos, ea_u, nea_u, ea_c, nea_c)],
      db_exp,
      by = "pos",
      all = FALSE,
      allow.cartesian = TRUE,
      sort = FALSE
    )
    if (nrow(cand) == 0) next

    cand[, rescue_rank := fifelse(ea_u == alt_allele & nea_u == ref, 1L,
                           fifelse(ea_u == ref & nea_u == alt_allele, 2L,
                           fifelse(ea_c == alt_allele & nea_c == ref, 3L,
                           fifelse(ea_c == ref & nea_c == alt_allele, 4L, 99L))))]

    cand <- cand[rescue_rank < 99L]
    if (nrow(cand) == 0) next

    best <- cand[order(row_id, rescue_rank)][, .SD[1], by = row_id]
    best <- best[, .(row_id, mapped_rsid_rescue = rsid, rescue_rank)]
    rescue_chunks[[length(rescue_chunks) + 1]] <- best
  }

  if (length(rescue_chunks) > 0) {
    rescue_dt <- rbindlist(rescue_chunks, use.names = TRUE, fill = TRUE)
    rescue_dt <- rescue_dt[!duplicated(row_id)]

    out <- merge(out, rescue_dt, by = "row_id", all.x = TRUE, sort = FALSE)

    out[!is.na(mapped_rsid_rescue) & (is.na(mapped_rsid) | mapped_rsid == ""), mapped_rsid := mapped_rsid_rescue]
    out[!is.na(mapped_rsid_rescue) & (is.na(map_status) | map_status == "unmapped"),
        map_status := fifelse(rescue_rank == 1L, "multiallelic_direct",
                     fifelse(rescue_rank == 2L, "multiallelic_swapped",
                     fifelse(rescue_rank == 3L, "multiallelic_strand_direct", "multiallelic_strand_swapped")))]

    rescued_count <- out[!is.na(mapped_rsid_rescue), .N]
    out[, c("mapped_rsid_rescue", "rescue_rank") := NULL]
  }
}

# Final mapped/unmapped after rescue
out[, snpid_fallback := fifelse(
  !is.na(chrom) & !is.na(pos) & !is.na(nea) & !is.na(ea),
  paste0(as.character(chrom), ":", as.character(pos), ":", as.character(nea), ":", as.character(ea)),
  NA_character_
)]

out[, snpid := fifelse(
  !is.na(mapped_rsid) & mapped_rsid != "",
  mapped_rsid,
  fifelse(!is.na(snpid_fallback) & snpid_fallback != "", snpid_fallback, snpid)
)]
out[, snpid_fallback := NULL]

mapped_count <- out[!is.na(mapped_rsid) & mapped_rsid != "", .N]
unmapped_dt <- out[is.na(mapped_rsid) | mapped_rsid == ""]
unmapped_count <- nrow(unmapped_dt)
cat(sprintf("Mapped rsID count: %d / %d\n", mapped_count, nrow(out)))
cat(sprintf("Multiallelic rescue mapped: %d\n", rescued_count))
cat("Mapping status summary:\n")
status_summary <- out[, .N, by = map_status][order(-N)]
print(status_summary)

# Limited fallback: for only a few unmapped variants, search by chr+pos in reference chr VCF.
lookup_n <- min(unmapped_pos_limit, unmapped_count)
pos_lookup_dt <- data.table()

if (lookup_n > 0) {
  lookup_subset <- copy(unmapped_dt[1:lookup_n])
  if (!"chrom_norm" %in% names(lookup_subset)) {
    lookup_subset[, chrom_norm := normalize_chr(chrom)]
  }

  lookup_chunks <- list()
  chr_for_lookup <- unique(lookup_subset[!is.na(chrom_norm) & !is.na(pos), chrom_norm])

  for (chr_norm in chr_for_lookup) {
    chr_rows <- lookup_subset[chrom_norm == chr_norm & !is.na(pos)]
    if (nrow(chr_rows) == 0) next

    vcf_file <- find_chr_vcf(opt$snpdb_root, build_info, chr_norm)
    if (is.na(vcf_file)) next

    pos_dt <- unique(chr_rows[, .(pos = as.integer(pos))])
    pos_dt <- pos_dt[!is.na(pos)]
    if (nrow(pos_dt) == 0) next
    pos_dt[, end := pos]

    reg_file <- file.path(opt$temp_dir, sprintf("unmapped_lookup_%s.tsv", chr_norm))
    db <- query_db_for_positions(vcf_file, chr_norm, pos_dt, reg_file)
    if (nrow(db) == 0) next

    db <- db[!is.na(rsid) & rsid != ".", .(pos = as.integer(pos), pos_only_rsid = rsid)]
    if (nrow(db) == 0) next

    db <- db[, .SD[1], by = pos]

    setkey(chr_rows, pos)
    setkey(db, pos)
    cand <- db[chr_rows]
    cand <- cand[, .(row_id, chrom, pos, ea, nea, map_status, pos_only_rsid)]
    lookup_chunks[[length(lookup_chunks) + 1]] <- cand
  }

  if (length(lookup_chunks) > 0) {
    pos_lookup_dt <- rbindlist(lookup_chunks, use.names = TRUE, fill = TRUE)
    pos_lookup_dt <- unique(pos_lookup_dt, by = "row_id")

    out <- merge(
      out,
      pos_lookup_dt[, .(row_id, pos_only_rsid)],
      by = "row_id",
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    out[, pos_only_rsid := NA_character_]
  }

  cat(sprintf("Unmapped chr+pos fallback lookup checked: %d variants\n", lookup_n))
  cat(sprintf("Unmapped chr+pos fallback hits: %d\n", out[is.na(mapped_rsid) | mapped_rsid == "", sum(!is.na(pos_only_rsid) & pos_only_rsid != "")]))
} else {
  out[, pos_only_rsid := NA_character_]
}

unmapped_dt <- out[is.na(mapped_rsid) | mapped_rsid == ""]
fwrite(unmapped_dt, unmapped_file, sep = "\t", na = "NA")

unmapped_with_pos_lookup <- out[(is.na(mapped_rsid) | mapped_rsid == "") & !is.na(pos_only_rsid) & pos_only_rsid != "",
                                .(row_id, chrom, pos, ea, nea, map_status, pos_only_rsid)]
if (lookup_n > 0) {
  fwrite(unmapped_with_pos_lookup, unmapped_pos_lookup_file, sep = "\t", na = "NA")
} else {
  fwrite(unmapped_with_pos_lookup[0], unmapped_pos_lookup_file, sep = "\t", na = "NA")
}

final_unmapped_count <- nrow(unmapped_dt)

# Remove helper columns not needed in final output
out[, c("row_id", "chrom_norm") := NULL]
fwrite(out, opt$output, sep = "\t", na = "NA")

summary_dt <- data.table(
  dataset = opt$dataset,
  build = opt$build,
  total_variants = nrow(out),
  mapped_rsid_variants = mapped_count,
  existing_snpid_variants = NA_integer_,
  mode = "mapped_with_snpdb"
)
fwrite(summary_dt, opt$summary, sep = "\t")

status_out <- copy(status_summary)
status_out[, `:=`(dataset = opt$dataset, build = opt$build)]
setcolorder(status_out, c("dataset", "build", "map_status", "N"))
fwrite(status_out, file.path(dirname(opt$summary), sprintf("%s_mapping_status.tsv", opt$dataset)), sep = "\t")

with(file(opt$done, "w"), {
  writeLines(c(
    sprintf("SNPdb mapping completed for %s", opt$dataset),
    sprintf("Build: %s", opt$build),
    sprintf("Input: %s", opt$input),
    sprintf("Output: %s", opt$output),
    sprintf("Summary: %s", opt$summary),
    sprintf("Unmapped list: %s", unmapped_file),
    sprintf("Unmapped chr+pos lookup list: %s", unmapped_pos_lookup_file),
    sprintf("Temp dir: %s", opt$temp_dir),
    sprintf("Total variants: %d", nrow(out)),
    sprintf("Mapped rsID variants: %d", mapped_count),
    sprintf("Mapped by multiallelic rescue: %d", rescued_count),
    sprintf("Unmapped variants: %d", final_unmapped_count),
    sprintf("Unmapped chr+pos fallback checked: %d", lookup_n),
    sprintf("Unmapped chr+pos fallback hits: %d", nrow(unmapped_with_pos_lookup))
  ))
})

cat(sprintf("Unmapped list written: %s (%d variants)\n", unmapped_file, final_unmapped_count))
cat(sprintf("Unmapped chr+pos lookup written: %s (%d candidate hits)\n", unmapped_pos_lookup_file, nrow(unmapped_with_pos_lookup)))

cat("Done.\n")
