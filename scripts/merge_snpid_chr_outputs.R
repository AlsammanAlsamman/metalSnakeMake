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
  make_option("--chrom-order", type = "character", default = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,X,Y"),
  make_option("--inputs-csv", type = "character"),
  make_option("--output", type = "character"),
  make_option("--done", type = "character"),
  make_option("--summary", type = "character"),
  make_option("--log", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("dataset", "output", "done", "summary", "log")
missing_required <- required[sapply(required, function(k) is.null(opt[[k]]) || !nzchar(opt[[k]]))]
if (length(missing_required) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_required, collapse = ", ")))
}

if (is.null(opt$`inputs-csv`) || !nzchar(opt$`inputs-csv`)) {
  stop("No per-chromosome inputs were provided to merge.")
}

dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$done), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$log), recursive = TRUE, showWarnings = FALSE)

chrom_order <- strsplit(opt$`chrom-order`, ",", fixed = TRUE)[[1]]
chrom_order <- trimws(chrom_order)

extract_chr <- function(path) {
  base <- basename(path)
  sub("^chr_([^.]+)\\..*$", "\\1", base)
}

input_paths <- strsplit(opt$`inputs-csv`, ",", fixed = TRUE)[[1]]
input_paths <- trimws(input_paths)
input_paths <- input_paths[nzchar(input_paths)]

input_dt <- data.table(path = input_paths)
input_dt[, chrom := extract_chr(path)]
input_dt[, ord := match(chrom, chrom_order)]
input_dt[is.na(ord), ord := .Machine$integer.max]
setorder(input_dt, ord, chrom)

log_msg("Merging ", nrow(input_dt), " per-chromosome files for dataset ", opt$dataset)

parts <- lapply(input_dt$path, function(p) {
  if (!file.exists(p)) {
    stop(sprintf("Missing per-chrom mapped file: %s", p))
  }
  fread(p, sep = "\t", showProgress = FALSE)
})

out <- rbindlist(parts, use.names = TRUE, fill = TRUE)
if ("row_id" %in% names(out)) {
  setorder(out, row_id)
}

mapped_count <- if ("mapped_rsid" %in% names(out)) out[!is.na(mapped_rsid) & mapped_rsid != "", .N] else 0L
existing_snpid_count <- if ("map_status" %in% names(out)) out[map_status == "existing_snpid", .N] else 0L

status_summary <- if ("map_status" %in% names(out)) {
  out[, .N, by = map_status][order(-N)]
} else {
  data.table(map_status = character(), N = integer())
}

unmapped_file <- file.path(dirname(opt$summary), sprintf("%s_unmapped.tsv", opt$dataset))
unmapped_pos_lookup_file <- file.path(dirname(opt$summary), sprintf("%s_unmapped_pos_lookup.tsv", opt$dataset))

if ("mapped_rsid" %in% names(out)) {
  unmapped_dt <- out[is.na(mapped_rsid) | mapped_rsid == ""]
} else {
  unmapped_dt <- out[0]
}

fwrite(unmapped_dt, unmapped_file, sep = "\t", na = "NA", quote = FALSE)

if (all(c("row_id", "chrom", "pos", "ea", "nea", "map_status", "pos_only_rsid") %in% names(out))) {
  unmapped_with_pos_lookup <- out[(is.na(mapped_rsid) | mapped_rsid == "") & !is.na(pos_only_rsid) & pos_only_rsid != "",
                                  .(row_id, chrom, pos, ea, nea, map_status, pos_only_rsid)]
} else {
  unmapped_with_pos_lookup <- data.table(row_id = integer(), chrom = character(), pos = integer(), ea = character(), nea = character(), map_status = character(), pos_only_rsid = character())
}
fwrite(unmapped_with_pos_lookup, unmapped_pos_lookup_file, sep = "\t", na = "NA", quote = FALSE)

out_write <- copy(out)
helper_cols <- c("row_id", "chrom_norm", "mapped_rsid", "map_status", "pos_only_rsid", "snpid_orig")
drop_cols <- intersect(helper_cols, names(out_write))
if (length(drop_cols) > 0) {
  out_write[, (drop_cols) := NULL]
}

join_artifact_cols <- grep("^i\\.", names(out_write), value = TRUE)
if (length(join_artifact_cols) > 0) {
  out_write[, (join_artifact_cols) := NULL]
}

strip_quotes <- function(x) {
  gsub('"', "", as.character(x), fixed = TRUE)
}

for (col in intersect(c("chrom", "ea", "nea", "snpid"), names(out_write))) {
  out_write[, (col) := strip_quotes(get(col))]
}

fwrite(out_write, opt$output, sep = "\t", na = "NA", quote = FALSE)

summary_dt <- data.table(
  dataset = opt$dataset,
  build = NA_character_,
  total_variants = nrow(out_write),
  mapped_rsid_variants = mapped_count,
  existing_snpid_variants = existing_snpid_count,
  mode = "mapped_with_snpdb_chr_jobs"
)
fwrite(summary_dt, opt$summary, sep = "\t", quote = FALSE)

status_out <- copy(status_summary)
if (nrow(status_out) > 0) {
  status_out[, `:=`(dataset = opt$dataset, build = NA_character_)]
  setcolorder(status_out, c("dataset", "build", "map_status", "N"))
  fwrite(status_out, file.path(dirname(opt$summary), sprintf("%s_mapping_status.tsv", opt$dataset)), sep = "\t", quote = FALSE)
} else {
  fwrite(data.table(dataset = character(), build = character(), map_status = character(), N = integer()), file.path(dirname(opt$summary), sprintf("%s_mapping_status.tsv", opt$dataset)), sep = "\t", quote = FALSE)
}

writeLines(c(
  sprintf("SNPdb mapping completed for %s", opt$dataset),
  sprintf("Output: %s", opt$output),
  sprintf("Summary: %s", opt$summary),
  sprintf("Unmapped list: %s", unmapped_file),
  sprintf("Unmapped chr+pos lookup list: %s", unmapped_pos_lookup_file),
  sprintf("Total variants: %d", nrow(out_write)),
  sprintf("Mapped rsID variants: %d", mapped_count),
  sprintf("Existing snpid variants: %d", existing_snpid_count),
  sprintf("Unmapped variants: %d", nrow(unmapped_dt)),
  sprintf("Unmapped chr+pos fallback hits: %d", nrow(unmapped_with_pos_lookup))
), con = opt$done)

log_msg("Merge complete. Total=", nrow(out_write), ", mapped=", mapped_count)
