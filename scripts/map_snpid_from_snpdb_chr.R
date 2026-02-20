#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

option_list <- list(
  make_option("--dataset", type = "character"),
  make_option("--chrom", type = "character"),
  make_option("--build", type = "character"),
  make_option("--input", type = "character"),
  make_option("--output", type = "character"),
  make_option("--done", type = "character"),
  make_option("--log", type = "character"),
  make_option("--snpdb-root", type = "character", default = "resources/SNPdb"),
  make_option("--mapping-dir", type = "character", default = "results/snpdmapped"),
  make_option("--bcftools", type = "character", default = "bcftools")
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

opt$snpdb_root <- get_opt(opt, c("snpdb_root", "snpdb-root", "snpdb.root"))
opt$mapping_dir <- get_opt(opt, c("mapping_dir", "mapping-dir", "mapping.dir"))

required <- c("dataset", "chrom", "build", "input", "output", "done", "log", "snpdb_root", "mapping_dir", "bcftools")
missing_required <- required[sapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]))]
if (length(missing_required) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_required, collapse = ", ")))
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$done), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$log), recursive = TRUE, showWarnings = FALSE)
dir.create(opt$mapping_dir, recursive = TRUE, showWarnings = FALSE)

normalize_chr <- function(x) {
  x <- toupper(as.character(x))
  x <- sub("^CHR", "", x)
  x[x == "M"] <- "MT"
  x
}

chrom_aliases <- function(chrom_norm) {
  ch <- toupper(as.character(chrom_norm))
  out <- c(ch)
  if (ch == "23") out <- c(out, "X")
  if (ch == "24") out <- c(out, "Y")
  if (ch == "25") out <- c(out, "MT", "M")
  if (ch == "X") out <- c(out, "23")
  if (ch == "Y") out <- c(out, "24")
  if (ch %in% c("MT", "M")) out <- c(out, "25")
  unique(out)
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
  stop(sprintf("Unsupported build '%s'.", build))
}

find_chr_vcf <- function(snpdb_root, build_info, chrom_norm) {
  aliases <- chrom_aliases(chrom_norm)
  chrom_candidates <- unique(c(aliases, paste0("chr", aliases), if (any(aliases %in% c("MT", "M"))) c("M", "chrM") else character(0)))
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

query_db_for_positions <- function(vcf_file, chrom_norm, pos_dt, regions_file, bcftools_cmd) {
  if (!exists(".vcf_fallback_cache", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".vcf_fallback_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  fallback_cache <- get(".vcf_fallback_cache", envir = .GlobalEnv, inherits = FALSE)

  bcftools_ok <- FALSE
  if (!is.null(bcftools_cmd) && nzchar(bcftools_cmd)) {
    bcftools_ok <- suppressWarnings(tryCatch(system2(bcftools_cmd, "--version", stdout = FALSE, stderr = FALSE) == 0, error = function(e) FALSE))
  }

  if (!bcftools_ok) {
    key <- gsub("[^A-Za-z0-9]", "_", normalizePath(vcf_file, winslash = "/", mustWork = FALSE))
    if (!exists(key, envir = fallback_cache, inherits = FALSE)) {
      log_msg("bcftools unavailable; using fallback full-scan for ", basename(vcf_file))
      pos_file <- tempfile(pattern = "snpdb_pos_", fileext = ".txt")
      out_file <- tempfile(pattern = "snpdb_scan_", fileext = ".tsv")
      on.exit(unlink(c(pos_file, out_file), force = TRUE), add = TRUE)

      fwrite(unique(pos_dt[, .(pos)]), pos_file, sep = "\t", col.names = FALSE)
      awk_script <- "BEGIN{FS=OFS=\"\\t\"} NR==FNR{p[$1]=1;next} /^#/{next} (($2+0) in p){print $1,$2,$3,$4,$5}"
      cmd <- sprintf("gzip -dc %s | awk '%s' %s - > %s", shQuote(vcf_file), awk_script, shQuote(pos_file), shQuote(out_file))
      status <- suppressWarnings(system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = FALSE))
      if (!is.numeric(status) || status != 0 || !file.exists(out_file) || file.info(out_file)$size == 0) {
        empty <- data.table(db_chrom = character(), pos = integer(), rsid = character(), ref = character(), alt = character())
        assign(key, empty, envir = fallback_cache)
      } else {
        db_all <- fread(out_file, sep = "\t", header = FALSE, showProgress = FALSE)
        if (nrow(db_all) == 0) {
          db_all <- data.table(db_chrom = character(), pos = integer(), rsid = character(), ref = character(), alt = character())
        } else {
          setnames(db_all, c("db_chrom", "pos", "rsid", "ref", "alt"))
          db_all[, pos := as.integer(pos)]
        }
        assign(key, db_all, envir = fallback_cache)
      }
    }

    db_cached <- get(key, envir = fallback_cache, inherits = FALSE)
    if (nrow(db_cached) == 0) {
      return(data.table(db_chrom = character(), pos = integer(), rsid = character(), ref = character(), alt = character()))
    }
    lookup_pos <- unique(as.integer(pos_dt$pos))
    lookup_pos <- lookup_pos[!is.na(lookup_pos)]
    if (length(lookup_pos) == 0) {
      return(data.table(db_chrom = character(), pos = integer(), rsid = character(), ref = character(), alt = character()))
    }
    return(db_cached[pos %in% lookup_pos])
  }

  aliases <- chrom_aliases(chrom_norm)
  contig_candidates <- unique(c(aliases, paste0("chr", aliases), if (any(aliases %in% c("MT", "M"))) c("M", "chrM") else character(0)))

  for (contig in contig_candidates) {
    reg <- copy(pos_dt)
    reg[, chrom := contig]
    setcolorder(reg, c("chrom", "pos", "end"))
    fwrite(reg, regions_file, sep = "\t", col.names = FALSE)

    cmd <- sprintf("%s query -R %s -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\n' %s", shQuote(bcftools_cmd), shQuote(regions_file), shQuote(vcf_file))
    db <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, showProgress = FALSE), error = function(e) data.table())
    if (nrow(db) > 0) {
      setnames(db, c("db_chrom", "pos", "rsid", "ref", "alt"))
      return(db)
    }
  }
  data.table(db_chrom = character(), pos = integer(), rsid = character(), ref = character(), alt = character())
}

run_multiallelic_rescue <- function(dt_chr, chrom_norm, vcf_file, mapping_dir, bcftools_cmd) {
  unmapped <- dt_chr[is.na(mapped_rsid) | mapped_rsid == ""]
  if (nrow(unmapped) == 0) return(dt_chr)

  pos_dt <- unique(unmapped[, .(pos = as.integer(pos))])
  pos_dt <- pos_dt[!is.na(pos)]
  if (nrow(pos_dt) == 0) return(dt_chr)
  pos_dt[, end := pos]

  reg_file <- file.path(mapping_dir, sprintf("chr_%s_multiallelic_rescue_regions.tsv", chrom_norm))
  db <- query_db_for_positions(vcf_file, chrom_norm, pos_dt, reg_file, bcftools_cmd)
  if (nrow(db) == 0) return(dt_chr)

  db <- db[!is.na(rsid) & rsid != ".", .(pos = as.integer(pos), rsid, ref = toupper(ref), alt = toupper(alt))]
  db <- db[!is.na(ref) & ref != "" & !is.na(alt) & alt != ""]
  if (nrow(db) == 0) return(dt_chr)

  db_exp <- db[, .(alt_allele = unlist(strsplit(alt, ",", fixed = TRUE))), by = .(pos, rsid, ref)]
  db_exp <- db_exp[!is.na(alt_allele) & alt_allele != "" & alt_allele != "."]
  if (nrow(db_exp) == 0) return(dt_chr)

  work <- copy(unmapped)
  work[, `:=`(ea_u = toupper(as.character(ea)), nea_u = toupper(as.character(nea)), ea_c = complement_allele(ea), nea_c = complement_allele(nea))]

  cand <- merge(work[, .(row_id, pos, ea_u, nea_u, ea_c, nea_c)], db_exp, by = "pos", all = FALSE, allow.cartesian = TRUE, sort = FALSE)
  if (nrow(cand) == 0) return(dt_chr)

  cand[, rescue_rank := fifelse(ea_u == alt_allele & nea_u == ref, 1L,
                         fifelse(ea_u == ref & nea_u == alt_allele, 2L,
                         fifelse(ea_c == alt_allele & nea_c == ref, 3L,
                         fifelse(ea_c == ref & nea_c == alt_allele, 4L, 99L))))]
  cand <- cand[rescue_rank < 99L]
  if (nrow(cand) == 0) return(dt_chr)

  best <- cand[order(row_id, rescue_rank)][, .SD[1], by = row_id][, .(row_id, mapped_rsid_rescue = rsid, rescue_rank)]
  out <- merge(dt_chr, best, by = "row_id", all.x = TRUE, sort = FALSE)
  out[!is.na(mapped_rsid_rescue) & (is.na(mapped_rsid) | mapped_rsid == ""), mapped_rsid := mapped_rsid_rescue]
  out[!is.na(mapped_rsid_rescue) & (is.na(map_status) | map_status == "unmapped"),
      map_status := fifelse(rescue_rank == 1L, "multiallelic_direct", fifelse(rescue_rank == 2L, "multiallelic_swapped", fifelse(rescue_rank == 3L, "multiallelic_strand_direct", "multiallelic_strand_swapped")))]
  out[, c("mapped_rsid_rescue", "rescue_rank") := NULL]
  out
}

log_msg("Dataset=", opt$dataset, " chrom=", opt$chrom)

chr_target <- normalize_chr(opt$chrom)
dt <- fread(opt$input, sep = "\t", showProgress = TRUE)

strip_quotes <- function(x) {
  gsub('"', "", as.character(x), fixed = TRUE)
}

for (col in intersect(c("chrom", "ea", "nea", "snpid"), names(dt))) {
  dt[, (col) := strip_quotes(get(col))]
}

required_cols <- c("chrom", "pos", "ea", "nea")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop(sprintf("Input missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

if (!"snpid" %in% names(dt)) dt[, snpid := NA_character_]

dt[, row_id := .I]
dt[, chrom_norm := normalize_chr(chrom)]
dt[, pos := suppressWarnings(as.integer(pos))]
dt_chr <- copy(dt[chrom_norm == chr_target])

if (nrow(dt_chr) == 0) {
  empty <- dt[0]
  empty[, `:=`(snpid_orig = character(), mapped_rsid = character(), map_status = character(), pos_only_rsid = character())]
  fwrite(empty, opt$output, sep = "\t", na = "NA", quote = FALSE)
  writeLines(c(sprintf("Dataset: %s", opt$dataset), sprintf("Chromosome: %s", chr_target), "Status: done (no variants)"), con = opt$done)
  quit(save = "no", status = 0)
}

dt_chr[, snpid_orig := snpid]
build_info <- build_to_keys(opt$build)
vcf_file <- find_chr_vcf(opt$snpdb_root, build_info, chr_target)

if (dt_chr[!is.na(snpid) & trimws(as.character(snpid)) != "", .N] == nrow(dt_chr)) {
  dt_chr[, `:=`(mapped_rsid = as.character(snpid), map_status = "existing_snpid", pos_only_rsid = NA_character_)]
} else if (is.na(vcf_file)) {
  dt_chr[, `:=`(mapped_rsid = NA_character_, map_status = "missing_snpdb_chr", pos_only_rsid = NA_character_)]
} else {
  pos_dt <- unique(dt_chr[, .(pos = as.integer(pos))])
  pos_dt <- pos_dt[!is.na(pos)]
  pos_dt[, end := pos]

  if (nrow(pos_dt) == 0) {
    dt_chr[, `:=`(mapped_rsid = NA_character_, map_status = "no_pos", pos_only_rsid = NA_character_)]
  } else {
    reg_file <- file.path(opt$mapping_dir, sprintf("chr_%s_regions.tsv", chr_target))
    db <- query_db_for_positions(vcf_file, chr_target, pos_dt, reg_file, opt$bcftools)

    if (nrow(db) == 0) {
      dt_chr[, `:=`(mapped_rsid = NA_character_, map_status = "no_db_match_pos", pos_only_rsid = NA_character_)]
    } else {
      db <- db[!is.na(rsid) & rsid != "."]
      if (nrow(db) == 0) {
        dt_chr[, `:=`(mapped_rsid = NA_character_, map_status = "no_valid_rsid", pos_only_rsid = NA_character_)]
      } else {
        db <- db[, .(pos = as.integer(pos), rsid, ref = toupper(ref), alt = toupper(alt))]
        db <- db[!grepl(",", alt)]

        dt_chr[, `:=`(ea_u = toupper(as.character(ea)), nea_u = toupper(as.character(nea)), ea_c = complement_allele(ea), nea_c = complement_allele(nea))]

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
        } else {
          best <- cands[order(row_id, rank)][, .SD[1], by = row_id][, .(row_id, mapped_rsid = rsid, map_rank = rank)]
          setkey(dt_chr, row_id)
          setkey(best, row_id)
          dt_chr <- best[dt_chr]
          dt_chr[, map_status := fifelse(is.na(mapped_rsid), "unmapped",
                                  fifelse(map_rank == 1L, "direct",
                                  fifelse(map_rank == 2L, "swapped",
                                  fifelse(map_rank == 3L, "strand_direct", "strand_swapped"))))]
          dt_chr[, map_rank := NULL]
        }

        dt_chr <- run_multiallelic_rescue(dt_chr, chr_target, vcf_file, opt$mapping_dir, opt$bcftools)
        dt_chr[, pos_only_rsid := NA_character_]

        unmapped <- dt_chr[is.na(mapped_rsid) | mapped_rsid == ""]
        if (nrow(unmapped) > 0) {
          pos_lookup <- unique(unmapped[!is.na(pos), .(pos = as.integer(pos))])
          pos_lookup <- pos_lookup[!is.na(pos)]
          if (nrow(pos_lookup) > 0) {
            pos_lookup[, end := pos]
            lookup_reg <- file.path(opt$mapping_dir, sprintf("chr_%s_unmapped_lookup_regions.tsv", chr_target))
            db_lookup <- query_db_for_positions(vcf_file, chr_target, pos_lookup, lookup_reg, opt$bcftools)
            if (nrow(db_lookup) > 0) {
              db_lookup <- db_lookup[!is.na(rsid) & rsid != ".", .(pos = as.integer(pos), pos_only_rsid = rsid)]
              if (nrow(db_lookup) > 0) {
                db_lookup <- db_lookup[, .SD[1], by = pos]
                setkey(unmapped, pos)
                setkey(db_lookup, pos)
                look <- db_lookup[unmapped][, .(row_id, pos_only_rsid)]
                look <- look[!is.na(row_id)]
                if (nrow(look) > 0) {
                  setkey(dt_chr, row_id)
                  setkey(look, row_id)
                  dt_chr <- look[dt_chr]
                }
              }
            }
          }
        }

        dt_chr[, c("ea_u", "nea_u", "ea_c", "nea_c") := NULL]
      }
    }
  }
}

dt_chr[, snpid_fallback := fifelse(!is.na(chrom) & !is.na(pos) & !is.na(nea) & !is.na(ea), paste0(as.character(chrom), ":", as.character(pos), ":", as.character(nea), ":", as.character(ea)), NA_character_)]

dt_chr[, snpid := fifelse(!is.na(mapped_rsid) & mapped_rsid != "", mapped_rsid, fifelse(!is.na(snpid_fallback) & snpid_fallback != "", snpid_fallback, snpid))]
dt_chr[, snpid_fallback := NULL]

for (col in intersect(c("chrom", "ea", "nea", "snpid"), names(dt_chr))) {
  dt_chr[, (col) := strip_quotes(get(col))]
}

if ("i.pos_only_rsid" %in% names(dt_chr)) {
  if (!"pos_only_rsid" %in% names(dt_chr)) {
    dt_chr[, pos_only_rsid := NA_character_]
  }
  dt_chr[is.na(pos_only_rsid) | pos_only_rsid == "", pos_only_rsid := as.character(`i.pos_only_rsid`)]
  dt_chr[, `i.pos_only_rsid` := NULL]
}

fwrite(dt_chr, opt$output, sep = "\t", na = "NA", quote = FALSE)

mapped_count <- dt_chr[!is.na(mapped_rsid) & mapped_rsid != "", .N]
rescued_count <- dt_chr[grepl("^multiallelic_", map_status), .N]
writeLines(c(
  sprintf("Dataset: %s", opt$dataset),
  sprintf("Chromosome: %s", chr_target),
  sprintf("Rows: %d", nrow(dt_chr)),
  sprintf("Mapped rsID rows: %d", mapped_count),
  sprintf("Multiallelic rescue rows: %d", rescued_count)
), con = opt$done)

log_msg("Chromosome ", chr_target, " done. rows=", nrow(dt_chr), ", mapped=", mapped_count)
