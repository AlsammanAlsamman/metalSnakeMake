#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list()

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (!startsWith(key, "--")) {
      stop(paste("Invalid argument:", key))
    }
    if (i == length(args)) {
      stop(paste("Missing value for", key))
    }
    out[[substring(key, 3)]] <- args[i + 1]
    i <- i + 2
  }

  required <- c(
    "combination", "datasets", "input_dir", "output_pdf",
    "output_summary", "output_table", "done", "log"
  )
  missing <- setdiff(required, names(out))
  if (length(missing) > 0) {
    stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
  }

  out
}

is_palindromic <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

log_msg <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}

read_dataset <- function(path, dataset_name) {
  if (!file.exists(path)) {
    stop(sprintf("Input file not found for dataset %s: %s", dataset_name, path))
  }

  required_cols <- c("chrom", "pos", "ea", "nea", "beta", "p")

  header_cols <- tryCatch(
    names(fread(path, sep = "\t", nrows = 0, showProgress = FALSE)),
    error = function(e) stop(sprintf("Failed reading header %s: %s", path, e$message))
  )

  missing_cols <- setdiff(required_cols, header_cols)
  if (length(missing_cols) > 0) {
    stop(sprintf("Dataset %s missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")))
  }

  df <- tryCatch(
    as.data.frame(
      fread(
        path,
        sep = "\t",
        select = required_cols,
        showProgress = FALSE
      )
    ),
    error = function(e) stop(sprintf("Failed reading %s: %s", path, e$message))
  )

  df <- df[!is.na(df$chrom) & !is.na(df$pos) & !is.na(df$ea) & !is.na(df$nea) & !is.na(df$beta) & !is.na(df$p), ]
  df$p <- suppressWarnings(as.numeric(df$p))
  df$beta <- suppressWarnings(as.numeric(df$beta))
  df$pos <- suppressWarnings(as.numeric(df$pos))
  df <- df[!is.na(df$p) & !is.na(df$beta) & !is.na(df$pos), ]
  df <- df[df$p < 5e-5, ]

  if (nrow(df) == 0) {
    return(df.frame())
  }

  df$key <- paste(df$chrom, df$pos, sep = ":")

  ord <- order(df$p)
  df <- df[ord, ]
  keep <- !duplicated(df$key)
  df <- df[keep, c("key", "chrom", "pos", "ea", "nea", "beta", "p")]

  rownames(df) <- NULL
  df
}

main <- function() {
  args <- parse_args()

  dir.create(dirname(args$output_pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$output_summary), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$output_table), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$done), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(args$log), recursive = TRUE, showWarnings = FALSE)

  datasets <- unlist(strsplit(args$datasets, ",", fixed = TRUE))
  datasets <- trimws(datasets)
  if (length(datasets) < 2) {
    stop("At least 2 datasets are required for correlation analysis")
  }

  log_msg(sprintf("Combination: %s", args$combination))
  log_msg(sprintf("Datasets: %s", paste(datasets, collapse = ", ")))

  data_map <- list()
  for (dataset_name in datasets) {
    input_file <- file.path(args$input_dir, paste0(dataset_name, ".tsv"))
    log_msg(sprintf("Reading dataset: %s (%s)", dataset_name, input_file))
    d <- read_dataset(input_file, dataset_name)
    log_msg(sprintf("Rows after p<5e-5 filter for %s: %d", dataset_name, nrow(d)))
    if (nrow(d) == 0) {
      stop(sprintf("No SNPs remaining for dataset %s after p<5e-5", dataset_name))
    }
    data_map[[dataset_name]] <- d
  }

  shared_keys <- Reduce(intersect, lapply(data_map, function(x) x$key))
  if (length(shared_keys) == 0) {
    stop("No shared SNPs found across all datasets at p<5e-5")
  }

  log_msg(sprintf("Shared SNP count before top-1000 selection: %d", length(shared_keys)))

  p_mat <- sapply(datasets, function(ds) {
    match_idx <- match(shared_keys, data_map[[ds]]$key)
    data_map[[ds]]$p[match_idx]
  })

  mean_p <- rowMeans(p_mat, na.rm = TRUE)
  ord <- order(mean_p)
  keep_n <- min(1000, length(shared_keys))
  selected_keys <- shared_keys[ord[seq_len(keep_n)]]
  log_msg(sprintf("Selected SNP count: %d", keep_n))

  reference_dataset <- datasets[1]
  ref_df <- data_map[[reference_dataset]]
  ref_df <- ref_df[match(selected_keys, ref_df$key), ]

  non_pal <- !is_palindromic(ref_df$ea, ref_df$nea)
  selected_keys <- selected_keys[non_pal]
  ref_df <- ref_df[non_pal, ]

  if (length(selected_keys) == 0) {
    stop("All selected shared SNPs are palindromic in reference dataset")
  }

  log_msg(sprintf("SNP count after removing palindromic reference SNPs: %d", length(selected_keys)))

  harmonized_long <- data.frame()

  for (dataset_name in datasets) {
    cur <- data_map[[dataset_name]]
    cur <- cur[match(selected_keys, cur$key), ]

    keep <- rep(TRUE, nrow(cur))
    beta_h <- cur$beta

    same <- toupper(cur$ea) == toupper(ref_df$ea) & toupper(cur$nea) == toupper(ref_df$nea)
    flip <- toupper(cur$ea) == toupper(ref_df$nea) & toupper(cur$nea) == toupper(ref_df$ea)
    valid <- same | flip

    keep[!valid] <- FALSE
    beta_h[flip] <- -beta_h[flip]

    tmp <- data.frame(
      key = cur$key,
      chrom = cur$chrom,
      pos = cur$pos,
      dataset = dataset_name,
      ea_ref = ref_df$ea,
      nea_ref = ref_df$nea,
      beta_h = beta_h,
      p = cur$p,
      keep = keep,
      stringsAsFactors = FALSE
    )
    harmonized_long <- rbind(harmonized_long, tmp)
  }

  keep_by_snp <- aggregate(keep ~ key, data = harmonized_long, FUN = function(x) all(x))
  valid_keys <- keep_by_snp$key[keep_by_snp$keep]

  harmonized_long <- harmonized_long[harmonized_long$key %in% valid_keys, ]
  harmonized_long$OR <- exp(harmonized_long$beta_h)

  if (length(valid_keys) < 10) {
    stop(sprintf("Too few harmonized shared SNPs after allele matching: %d", length(valid_keys)))
  }

  log_msg(sprintf("Final harmonized shared SNP count: %d", length(valid_keys)))

  ref_or <- harmonized_long[harmonized_long$dataset == reference_dataset, c("key", "OR")]
  colnames(ref_or)[2] <- "OR_ref"

  merged_all <- merge(harmonized_long, ref_or, by = "key", all.x = TRUE)
  plot_df <- merged_all[merged_all$dataset != reference_dataset, ]

  p <- ggplot(plot_df, aes(x = OR_ref, y = OR, color = dataset)) +
    geom_point(alpha = 0.45, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    labs(
      title = paste0("OR Concordance by Dataset: ", args$combination),
      subtitle = paste0("Shared SNPs (p<5e-5), top n=", length(valid_keys), ", reference=", reference_dataset),
      x = paste0("OR (reference: ", reference_dataset, ")"),
      y = "OR (dataset)",
      color = "Dataset"
    ) +
    theme_bw(base_size = 12)

  ggsave(args$output_pdf, p, width = 10, height = 7, units = "in")

  wide <- reshape(
    harmonized_long[, c("key", "dataset", "OR")],
    idvar = "key",
    timevar = "dataset",
    direction = "wide"
  )

  direction_mat <- as.data.frame(wide)
  or_cols <- grep("^OR\\.", names(direction_mat), value = TRUE)
  dir_bin <- as.data.frame(lapply(direction_mat[, or_cols, drop = FALSE], function(x) x > 1))
  majority_dir <- rowMeans(dir_bin, na.rm = TRUE) >= 0.5

  summary_rows <- list()
  for (dataset_name in datasets) {
    if (dataset_name == reference_dataset) {
      corr_val <- 1
      slope <- 1
    } else {
      ds_df <- merged_all[merged_all$dataset == dataset_name, c("OR_ref", "OR")]
      corr_val <- suppressWarnings(cor(ds_df$OR_ref, ds_df$OR, use = "complete.obs", method = "pearson"))
      lm_fit <- lm(OR ~ OR_ref, data = ds_df)
      slope <- unname(coef(lm_fit)[2])
    }

    ds_or_col <- paste0("OR.", dataset_name)
    ds_vec <- direction_mat[[ds_or_col]]
    disagree_majority <- mean((ds_vec > 1) != majority_dir, na.rm = TRUE)
    disagree_ref <- mean((ds_vec > 1) != (direction_mat[[paste0("OR.", reference_dataset)]] > 1), na.rm = TRUE)

    summary_rows[[dataset_name]] <- data.frame(
      combination = args$combination,
      reference_dataset = reference_dataset,
      dataset = dataset_name,
      n_snps_used = sum(!is.na(ds_vec)),
      pearson_with_ref = round(corr_val, 4),
      slope_with_ref = round(slope, 4),
      median_or = round(median(ds_vec, na.rm = TRUE), 4),
      disagree_with_ref_rate = round(disagree_ref, 4),
      disagree_with_majority_rate = round(disagree_majority, 4),
      flagged_different = ifelse((disagree_majority > 0.40) | (corr_val < 0.30), "YES", "NO"),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(-summary_df$disagree_with_majority_rate, summary_df$pearson_with_ref), ]

  write.table(summary_df, args$output_summary, sep = "\t", quote = FALSE, row.names = FALSE)

  table_out <- unique(harmonized_long[, c("key", "chrom", "pos", "ea_ref", "nea_ref")])
  for (dataset_name in datasets) {
    ds_out <- harmonized_long[harmonized_long$dataset == dataset_name, c("key", "beta_h", "OR", "p")]
    colnames(ds_out) <- c(
      "key",
      paste0("beta_h_", dataset_name),
      paste0("OR_", dataset_name),
      paste0("p_", dataset_name)
    )
    table_out <- merge(table_out, ds_out, by = "key", all.x = TRUE)
  }
  write.table(table_out, args$output_table, sep = "\t", quote = FALSE, row.names = FALSE)

  done_lines <- c(
    sprintf("Combination: %s", args$combination),
    sprintf("Reference dataset: %s", reference_dataset),
    sprintf("Final harmonized shared SNPs: %d", length(valid_keys)),
    sprintf("Output PDF: %s", args$output_pdf),
    sprintf("Output summary: %s", args$output_summary),
    sprintf("Output table: %s", args$output_table)
  )
  writeLines(done_lines, con = args$done)

  log_msg("Completed successfully")
}

tryCatch(
  main(),
  error = function(e) {
    cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
    quit(status = 1)
  }
)
