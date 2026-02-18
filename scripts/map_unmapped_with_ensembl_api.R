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
  make_option("--build", type = "character"),
  make_option("--mapped-input", type = "character"),
  make_option("--unmapped-input", type = "character"),
  make_option("--output", type = "character"),
  make_option("--done", type = "character"),
  make_option("--summary", type = "character"),
  make_option("--log", type = "character"),
  make_option("--rate-per-second", type = "double", default = 10)
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

opt$mapped_input <- get_opt(opt, c("mapped_input", "mapped-input", "mapped.input"))
opt$unmapped_input <- get_opt(opt, c("unmapped_input", "unmapped-input", "unmapped.input"))
opt$rate_per_second <- as.numeric(get_opt(opt, c("rate_per_second", "rate-per-second", "rate.per.second")))
if (is.na(opt$rate_per_second) || opt$rate_per_second <= 0) {
  opt$rate_per_second <- 10
}

required <- c("dataset", "build", "mapped_input", "unmapped_input", "output", "done", "summary", "log")
missing_required <- required[sapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]))]
if (length(missing_required) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_required, collapse = ", ")))
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$done), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$log), recursive = TRUE, showWarnings = FALSE)

strip_quotes <- function(x) {
  gsub('"', "", as.character(x), fixed = TRUE)
}

normalize_chr <- function(x) {
  y <- toupper(as.character(x))
  y <- sub("^CHR", "", y)
  y[y == "M"] <- "MT"
  y
}

api_base_url <- function(build) {
  b <- tolower(as.character(build))
  if (b %in% c("hg19", "grch37", "37")) {
    return("https://grch37.rest.ensembl.org")
  }
  if (b %in% c("hg38", "grch38", "38")) {
    return("https://rest.ensembl.org")
  }
  stop(sprintf("Unsupported build '%s'. Expected hg19/GRCh37 or hg38/GRCh38.", build))
}

parse_ensembl_yaml_html <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) {
    return(data.table(rsid = character(), alleles = character()))
  }

  pre <- txt
  if (grepl("<pre>", pre, fixed = TRUE)) {
    pre <- sub("(?s)^.*<pre>", "", pre, perl = TRUE)
    pre <- sub("(?s)</pre>.*$", "", pre, perl = TRUE)
  }

  lines <- unlist(strsplit(pre, "\n", fixed = TRUE))
  lines <- trimws(lines, which = "right")
  if (length(lines) == 0) {
    return(data.table(rsid = character(), alleles = character()))
  }

  records <- list()
  current_id <- NA_character_
  current_alleles <- character(0)
  in_alleles <- FALSE

  flush_record <- function() {
    if (!is.na(current_id) && nzchar(current_id)) {
      records[[length(records) + 1]] <<- data.table(
        rsid = current_id,
        alleles = paste(unique(current_alleles), collapse = ",")
      )
    }
  }

  for (ln in lines) {
    l <- trimws(ln)
    if (!nzchar(l) || l %in% c("---", "-")) next

    if (grepl("^alleles:\\s*$", l)) {
      in_alleles <- TRUE
      next
    }

    if (grepl("^id:\\s*", l)) {
      flush_record()
      current_id <- sub("^id:\\s*", "", l)
      current_id <- trimws(current_id)
      current_alleles <- character(0)
      in_alleles <- FALSE
      next
    }

    if (in_alleles && grepl("^-\\s*", l)) {
      a <- sub("^-\\s*", "", l)
      a <- toupper(trimws(a))
      if (nzchar(a)) {
        current_alleles <- c(current_alleles, a)
      }
      next
    }

    if (in_alleles && grepl(":", l, fixed = TRUE) && !grepl("^-\\s*", l)) {
      in_alleles <- FALSE
    }
  }

  flush_record()

  if (length(records) == 0) {
    return(data.table(rsid = character(), alleles = character()))
  }

  out <- rbindlist(records, use.names = TRUE, fill = TRUE)
  out <- out[!is.na(rsid) & nzchar(rsid)]
  out
}

query_ensembl_variation <- function(base_url, chrom, pos) {
  url <- sprintf("%s/overlap/region/human/%s:%d-%d?feature=variation", base_url, chrom, pos, pos)
  raw <- tryCatch(
    system2("curl", args = c("-s", url), stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )

  txt <- paste(raw, collapse = "\n")
  if (!nzchar(trimws(txt)) || grepl("No Data", txt, fixed = TRUE)) {
    return(data.table(rsid = character(), alleles = character()))
  }

  parse_ensembl_yaml_html(txt)
}

choose_rsid <- function(dt_api, ea, nea) {
  if (nrow(dt_api) == 0) {
    return(list(rsid = NA_character_, mapping = NA_character_))
  }

  ea_u <- toupper(as.character(ea))
  nea_u <- toupper(as.character(nea))

  dt_work <- copy(dt_api)
  dt_work[, allele_vec := strsplit(alleles, ",", fixed = TRUE)]
  dt_work[, allele_match := vapply(allele_vec, function(v) {
    vv <- toupper(trimws(v))
    any(vv %in% c(ea_u, nea_u))
  }, logical(1))]
  dt_work[, is_rsid := grepl("^rs", rsid, ignore.case = TRUE)]

  hit_allele <- dt_work[allele_match == TRUE]
  if (nrow(hit_allele) > 0) {
    hit_allele <- hit_allele[order(-is_rsid)]
    return(list(rsid = hit_allele$rsid[[1]], mapping = "allele"))
  }

  hit_pos <- dt_work[order(-is_rsid)]
  return(list(rsid = hit_pos$rsid[[1]], mapping = "position"))
}

log_msg("============================================================")
log_msg("API remap for unmapped SNPs")
log_msg("============================================================")
log_msg("Dataset: ", opt$dataset)
log_msg("Build: ", opt$build)
log_msg("Mapped input: ", opt$mapped_input)
log_msg("Unmapped input: ", opt$unmapped_input)
log_msg("Rate/sec: ", sprintf("%.2f", opt$rate_per_second))

if (!file.exists(opt$mapped_input)) stop(sprintf("Mapped input not found: %s", opt$mapped_input))
if (!file.exists(opt$unmapped_input)) stop(sprintf("Unmapped input not found: %s", opt$unmapped_input))

mapped <- fread(opt$mapped_input, sep = "\t", showProgress = TRUE)
unmapped <- fread(opt$unmapped_input, sep = "\t", showProgress = TRUE)

for (col in intersect(c("chrom", "ea", "nea", "snpid"), names(mapped))) {
  mapped[, (col) := strip_quotes(get(col))]
}
for (col in intersect(c("chrom", "ea", "nea", "snpid"), names(unmapped))) {
  unmapped[, (col) := strip_quotes(get(col))]
}

if (!all(c("chrom", "pos", "ea", "nea") %in% names(unmapped))) {
  stop("Unmapped input must include columns: chrom, pos, ea, nea")
}
if (!all(c("chrom", "pos", "ea", "nea", "snpid") %in% names(mapped))) {
  stop("Mapped input must include columns: chrom, pos, ea, nea, snpid")
}

unmapped[, chrom_norm := normalize_chr(chrom)]
unmapped[, pos := suppressWarnings(as.integer(pos))]
unmapped <- unmapped[!is.na(chrom_norm) & !is.na(pos)]

cands <- unique(unmapped[, .(chrom = as.character(chrom_norm), pos = as.integer(pos), ea = as.character(ea), nea = as.character(nea))])

if (nrow(cands) == 0) {
  mapped[, mapping := NA_character_]
  fwrite(mapped, opt$output, sep = "\t", na = "NA", quote = FALSE)
  summary_dt <- data.table(
    dataset = opt$dataset,
    build = opt$build,
    total_unmapped_candidates = 0L,
    api_mapped_total = 0L,
    api_mapped_allele = 0L,
    api_mapped_position = 0L,
    api_no_result = 0L,
    output_rows = nrow(mapped)
  )
  fwrite(summary_dt, opt$summary, sep = "\t", quote = FALSE)
  writeLines(c(
    sprintf("API remap completed for %s", opt$dataset),
    "No unmapped candidates found.",
    sprintf("Output: %s", opt$output)
  ), con = opt$done)
  quit(save = "no", status = 0)
}

base_url <- api_base_url(opt$build)
sleep_sec <- 1 / opt$rate_per_second

map_list <- vector("list", nrow(cands))
for (i in seq_len(nrow(cands))) {
  row <- cands[i]

  dt_api <- query_ensembl_variation(base_url, row$chrom[[1]], row$pos[[1]])
  picked <- choose_rsid(dt_api, row$ea[[1]], row$nea[[1]])

  map_list[[i]] <- data.table(
    chrom = row$chrom[[1]],
    pos = row$pos[[1]],
    ea = row$ea[[1]],
    nea = row$nea[[1]],
    api_rsid = picked$rsid,
    api_mapping = picked$mapping
  )

  if (i < nrow(cands)) Sys.sleep(sleep_sec)
}

api_map <- rbindlist(map_list, use.names = TRUE, fill = TRUE)
api_map <- api_map[!is.na(api_rsid) & nzchar(api_rsid)]

mapped[, chrom_norm := normalize_chr(chrom)]
mapped[, pos := suppressWarnings(as.integer(pos))]

if (nrow(api_map) > 0) {
  setkey(mapped, chrom_norm, pos, ea, nea)
  setkey(api_map, chrom, pos, ea, nea)
  joined <- api_map[mapped]

  if (!"mapping" %in% names(joined)) {
    joined[, mapping := NA_character_]
  }

  joined[!is.na(api_rsid) & nzchar(api_rsid), `:=`(
    snpid = api_rsid,
    mapping = fifelse(!is.na(api_mapping), api_mapping, "position")
  )]

  joined[, c("api_rsid", "api_mapping", "chrom") := NULL]
  setnames(joined, "chrom_norm", "chrom")
  mapped_out <- joined
} else {
  mapped[, mapping := NA_character_]
  mapped_out <- mapped
}

for (col in intersect(c("chrom", "ea", "nea", "snpid"), names(mapped_out))) {
  mapped_out[, (col) := strip_quotes(get(col))]
}

mapped_out[, intersect(c("chrom_norm"), names(mapped_out)) := NULL]
fwrite(mapped_out, opt$output, sep = "\t", na = "NA", quote = FALSE)

api_mapped_total <- if ("mapping" %in% names(mapped_out)) mapped_out[!is.na(mapping) & nzchar(mapping), .N] else 0L
api_mapped_allele <- if ("mapping" %in% names(mapped_out)) mapped_out[mapping == "allele", .N] else 0L
api_mapped_position <- if ("mapping" %in% names(mapped_out)) mapped_out[mapping == "position", .N] else 0L
api_no_result <- nrow(cands) - nrow(api_map)
if (api_no_result < 0) api_no_result <- 0L

summary_dt <- data.table(
  dataset = opt$dataset,
  build = opt$build,
  total_unmapped_candidates = nrow(cands),
  api_mapped_total = api_mapped_total,
  api_mapped_allele = api_mapped_allele,
  api_mapped_position = api_mapped_position,
  api_no_result = api_no_result,
  output_rows = nrow(mapped_out)
)
fwrite(summary_dt, opt$summary, sep = "\t", quote = FALSE)

writeLines(c(
  sprintf("API remap completed for %s", opt$dataset),
  sprintf("Build: %s", opt$build),
  sprintf("Candidates queried: %d", nrow(cands)),
  sprintf("API mapped total: %d", api_mapped_total),
  sprintf(" - allele match: %d", api_mapped_allele),
  sprintf(" - position only: %d", api_mapped_position),
  sprintf("No result from API: %d", api_no_result),
  sprintf("Output: %s", opt$output),
  sprintf("Summary: %s", opt$summary)
), con = opt$done)

log_msg("Done.")
