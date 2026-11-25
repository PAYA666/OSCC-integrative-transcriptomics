suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
first_perm <- if (length(args) >= 1L) as.integer(args[1]) else 1L
last_perm <- if (length(args) >= 2L) as.integer(args[2]) else 100L
stopifnot(1L <= first_perm, first_perm <= last_perm, last_perm <= 1000L)

meta <- fread("../01_baseline/results/sample_manifest.tsv", na.strings = "NA")[remove_base_and_mixed==TRUE]
v <- meta[sample_type == "Primary Tumor" & !is.na(OS_event),
          .(patient_barcode, sample_barcode, OS_time, OS_event)]
stopifnot(nrow(v) == 211L, !anyDuplicated(v$patient_barcode))

raw <- readRDS(file.path(Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA")),
                         "TCGA_OSCC_Nature8_raw_count_matrix.rds"))
raw <- raw[raw$gene_type == "protein_coding", ]
cnt <- as.matrix(raw[, meta$sample_barcode])
rownames(cnt) <- raw$gene_id
sym <- setNames(raw$gene_name, raw$gene_id)
storage.mode(cnt) <- "integer"
rm(raw)

trans <- function(q, dispersion) {
  a <- dispersion[["asymptDisp"]]
  e <- dispersion[["extraPois"]]
  log2((1 + e + 2 * a * q +
          2 * sqrt(a * q * (1 + e + a * q))) / (4 * a))
}

fold_files <- list.files("checkpoints", pattern = "^r[0-9]{2}f[0-9]{2}\\.rds$",
                         full.names = TRUE)
stopifnot(length(fold_files) == 100L)
folds <- lapply(fold_files, readRDS)
stopifnot(all(vapply(folds, function(z) length(z$candidate) > 1L, logical(1))))

prepared <- lapply(folds, function(a) {
  patients <- a$inner$patient
  train_meta <- v[match(patients, patient_barcode)]
  stopifnot(nrow(train_meta) == nrow(a$inner),
            !anyNA(train_meta$patient_barcode))
  genes <- intersect(a$candidate, rownames(cnt))
  samples <- train_meta$sample_barcode
  sf <- a$size_factors[samples]
  stopifnot(length(sf) == length(samples), all(is.finite(sf)),
            all(genes %in% names(a$screen_p)))
  q <- sweep(cnt[genes, samples, drop = FALSE], 2L, sf, "/")
  x <- t(trans(q, a$dispersion))
  colnames(x) <- genes
  stopifnot(nrow(x) == nrow(a$inner), all(is.finite(x)))
  list(x = x, time = train_meta$OS_time, event = train_meta$OS_event,
       patient = patients, original_inner = a$inner$inner_fold,
       candidate = genes, symbols = sym[genes], seed = a$seed,
       repeat_id = a$repeat_id, fold_id = a$fold_id)
})

cox_score_p <- function(x, time, event) {
  ord <- order(time)
  n <- length(time)
  reverse_x <- x[ord[n:1L], , drop = FALSE]
  reverse_x2 <- reverse_x * reverse_x
  cum_x <- matrixStats::colCumsums(reverse_x)
  cum_x2 <- matrixStats::colCumsums(reverse_x2)
  event_times <- sort(unique(time[event == 1L]))
  starts <- match(event_times, time[ord])
  ridx <- n - starts + 1L
  risk_n <- n - starts + 1L
  risk_x <- cum_x[ridx, , drop = FALSE]
  risk_x2 <- cum_x2[ridx, , drop = FALSE]
  event_rows <- which(event == 1L)
  group <- match(time[event_rows], event_times)
  d <- tabulate(group, nbins = length(event_times))
  event_x_raw <- rowsum(x[event_rows, , drop = FALSE], group, reorder = FALSE)
  event_x2_raw <- rowsum(x[event_rows, , drop = FALSE]^2, group, reorder = FALSE)
  event_order <- match(seq_along(d), as.integer(rownames(event_x_raw)))
  stopifnot(!anyNA(event_order))
  event_x <- event_x_raw[event_order, , drop = FALSE]
  event_x2 <- event_x2_raw[event_order, , drop = FALSE]
  mean_x <- risk_x / risk_n
  mean_x2 <- risk_x2 / risk_n
  score <- colSums(event_x) - colSums(mean_x * d)
  variance <- colSums((mean_x2 - mean_x^2) * d)

  tied <- which(d > 1L)
  for (i in tied) {
    old_mean <- mean_x[i, ]
    old_variance <- d[i] * (mean_x2[i, ] - old_mean^2)
    efron_mean_sum <- numeric(ncol(x))
    efron_variance_sum <- numeric(ncol(x))
    for (l in 0:(d[i] - 1L)) {
      frac <- l / d[i]
      denom <- risk_n[i] - l
      m <- (risk_x[i, ] - frac * event_x[i, ]) / denom
      m2 <- (risk_x2[i, ] - frac * event_x2[i, ]) / denom
      efron_mean_sum <- efron_mean_sum + m
      efron_variance_sum <- efron_variance_sum + pmax(0, m2 - m^2)
    }
    score <- score + d[i] * old_mean - efron_mean_sum
    variance <- variance - old_variance + efron_variance_sum
  }
  stat <- score^2 / variance
  pchisq(stat, df = 1L, lower.tail = FALSE)
}

screen_and_select <- function(d, time, event, perm_id) {
  p <- cox_score_p(d$x, time, event)
  names(p) <- colnames(d$x)
  use <- names(p)[which(is.finite(p) & p < 0.05)]
  if (length(use) < 2L) return(character())
  inner <- integer(length(event))
  set.seed(d$seed + 100000L * perm_id)
  for (e in 0:1) {
    idx <- which(event == e)
    inner[idx] <- sample(rep_len(1:10, length(idx)))
  }
  cv <- cv.glmnet(d$x[, use, drop = FALSE], Surv(time, event),
                  family = "cox", alpha = 1, standardize = TRUE,
                  foldid = inner, type.measure = "deviance")
  cf <- as.matrix(coef(cv, s = "lambda.min"))[, 1L]
  names(cf)[is.finite(cf) & cf != 0]
}

screen_comparison <- rbindlist(lapply(seq_along(prepared), function(i) {
  d <- prepared[[i]]
  p_score <- cox_score_p(d$x, d$time, d$event)
  p_wald <- folds[[i]]$screen_p[names(p_score)]
  data.table(repeat_id = d$repeat_id, fold_id = d$fold_id,
             n_candidates = length(p_score),
             wald_selected = sum(p_wald < 0.05, na.rm = TRUE),
             score_selected = sum(p_score < 0.05, na.rm = TRUE),
             threshold_agreement = mean((p_wald < 0.05) == (p_score < 0.05),
                                        na.rm = TRUE),
             selected_overlap = sum(p_wald < 0.05 & p_score < 0.05,
                                    na.rm = TRUE))
}))
fwrite(screen_comparison, "results/selection_screen_comparison.tsv", sep = "\t")
cat("Wald/score screen threshold agreement:",
    round(mean(screen_comparison$threshold_agreement), 4), "\n")
flush.console()

out <- "results/site_selection_null_permutations"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
final8 <- c("CTSG", "P4HA1", "NKAIN1", "LRRC38", "SPSB4", "SPINK6",
            "TNFRSF4", "DEFB1")

for (perm_id in seq.int(first_perm, last_perm)) {
  path <- file.path(out, sprintf("perm_%04d.tsv", perm_id))
  if (file.exists(path)) next
  set.seed(20261201L + perm_id)
  shuffled <- sample.int(nrow(v))
  perm_time <- v$OS_time[shuffled]
  perm_event <- v$OS_event[shuffled]
  selected <- lapply(prepared, function(d) {
    ii <- match(d$patient, v$patient_barcode)
    stopifnot(!anyNA(ii))
    screen_and_select(d, perm_time[ii], perm_event[ii], perm_id)
  })
  counts <- table(factor(unlist(lapply(seq_along(selected), function(i) {
    sym[selected[[i]]]
  })), levels = final8))
  result <- data.table(permutation = perm_id, gene = final8,
                       selected_outer_models = as.integer(counts),
                       selection_frequency = as.integer(counts) / length(folds))
  tmp <- paste0(path, ".tmp")
  fwrite(result, tmp, sep = "\t")
  stopifnot(file.rename(tmp, path))
  cat("Completed permutation", perm_id, "of", last_perm,
      "; any Final8 gene >=0.80:",
      any(result$selection_frequency >= 0.80), "\n")
  flush.console()
}

files <- list.files(out, pattern = "^perm_[0-9]{4}\\.tsv$", full.names = TRUE)
all <- rbindlist(lapply(files, fread))
summary <- all[, .(
  permutations = uniqueN(permutation),
  median_selection_frequency = median(selection_frequency),
  q95_selection_frequency = unname(quantile(selection_frequency, 0.95)),
  permutations_ge_0_80 = sum(selection_frequency >= 0.80)
), by = gene]
any_gene <- all[, .(any_final8_gene_ge_0_80 = any(selection_frequency >= 0.80)),
                by = permutation]
overall <- data.table(
  permutations = nrow(any_gene),
  permutations_with_any_final8_gene_ge_0_80 = sum(any_gene$any_final8_gene_ge_0_80),
  null_probability_any_final8_gene_ge_0_80 = mean(any_gene$any_final8_gene_ge_0_80)
)
  fwrite(summary, "results/site_selection_null_summary.tsv", sep = "\t")
  fwrite(overall, "results/site_selection_null_overall.tsv", sep = "\t")
capture.output(sessionInfo(), file = "results/site_selection_null_session.txt")
