suppressPackageStartupMessages({
  library(DESeq2)
  library(glmnet)
  library(survival)
  library(data.table)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
n_subsamples <- if (length(args) >= 1L) as.integer(args[1]) else 100L
stopifnot(n_subsamples >= 1L, n_subsamples <= 1000L)

meta <- fread("../01_baseline/results/sample_manifest.tsv", na.strings = "NA")
survival_cohort <- meta[sample_type == "Primary Tumor" & !is.na(OS_event),
                        .(patient_barcode, sample_barcode, OS_time, OS_event)]
stopifnot(nrow(survival_cohort) == 297L,
          !anyDuplicated(survival_cohort$patient_barcode))

frozen_pool <- fread(paste0(Sys.getenv("OSCC_DATA"), "/TCGA_OSCC_Nature8_Stage6F_DEG_nonGrey_candidates.tsv"))
candidate_ids <- frozen_pool$gene_id
stopifnot(length(candidate_ids) == 1543L, !anyDuplicated(candidate_ids))

raw <- readRDS(file.path(Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA")),
                         "TCGA_OSCC_Nature8_raw_count_matrix.rds"))
raw <- raw[raw$gene_type == "protein_coding", ]
cnt <- as.matrix(raw[, meta$sample_barcode])
rownames(cnt) <- raw$gene_id
sym <- setNames(raw$gene_name, raw$gene_id)
storage.mode(cnt) <- "integer"
rm(raw)
keep <- rowSums(cnt >= 10L) >= 10L
cnt <- cnt[keep, , drop = FALSE]
typ <- factor(ifelse(meta$sample_type == "Primary Tumor", "Tumor", "Normal"),
              levels = c("Normal", "Tumor"))
dds <- DESeqDataSetFromMatrix(cnt,
  data.frame(type = typ, row.names = colnames(cnt)), design = ~ type)
dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
dispersion <- attr(dispersionFunction(dds), "coefficients")
norm <- counts(dds, normalized = TRUE)
a <- dispersion[["asymptDisp"]]
e <- dispersion[["extraPois"]]
vst <- log2((1 + e + 2 * a * norm +
               2 * sqrt(a * norm * (1 + e + a * norm))) / (4 * a))
rm(dds, norm, cnt)

candidate_ids <- intersect(candidate_ids, rownames(vst))
stopifnot(length(candidate_ids) >= 1500L)
x <- t(vst[candidate_ids, survival_cohort$sample_barcode, drop = FALSE])
rownames(x) <- survival_cohort$patient_barcode
stopifnot(nrow(x) == 297L, all(is.finite(x)))

cox_score_p <- function(x, time, event) {
  ord <- order(time)
  n <- length(time)
  reverse_x <- x[ord[n:1L], , drop = FALSE]
  cum_x <- matrixStats::colCumsums(reverse_x)
  cum_x2 <- matrixStats::colCumsums(reverse_x * reverse_x)
  event_times <- sort(unique(time[event == 1L]))
  starts <- match(event_times, time[ord])
  risk_n <- n - starts + 1L
  ridx <- risk_n
  risk_x <- cum_x[ridx, , drop = FALSE]
  risk_x2 <- cum_x2[ridx, , drop = FALSE]
  event_rows <- which(event == 1L)
  group <- match(time[event_rows], event_times)
  d <- tabulate(group, nbins = length(event_times))
  event_x <- rowsum(x[event_rows, , drop = FALSE], group, reorder = FALSE)
  event_x2 <- rowsum(x[event_rows, , drop = FALSE]^2, group, reorder = FALSE)
  mean_x <- risk_x / risk_n
  mean_x2 <- risk_x2 / risk_n
  score <- colSums(event_x) - colSums(mean_x * d)
  variance <- colSums((mean_x2 - mean_x^2) * d)
  for (i in which(d > 1L)) {
    old_mean <- mean_x[i, ]
    old_variance <- d[i] * (mean_x2[i, ] - old_mean^2)
    mean_sum <- numeric(ncol(x))
    variance_sum <- numeric(ncol(x))
    for (l in 0:(d[i] - 1L)) {
      frac <- l / d[i]
      denom <- risk_n[i] - l
      m <- (risk_x[i, ] - frac * event_x[i, ]) / denom
      m2 <- (risk_x2[i, ] - frac * event_x2[i, ]) / denom
      mean_sum <- mean_sum + m
      variance_sum <- variance_sum + pmax(0, m2 - m^2)
    }
    score <- score + d[i] * old_mean - mean_sum
    variance <- variance - old_variance + variance_sum
  }
  pchisq(score^2 / variance, df = 1L, lower.tail = FALSE)
}

final8 <- c("CTSG", "P4HA1", "NKAIN1", "LRRC38", "SPSB4", "SPINK6",
            "TNFRSF4", "DEFB1")
set.seed(20261121L)
records <- vector("list", n_subsamples)
for (b in seq_len(n_subsamples)) {
  set.seed(20261121L + b)
  idx <- sample.int(nrow(survival_cohort), floor(nrow(survival_cohort) / 2L),
                    replace = FALSE)
  ytime <- survival_cohort$OS_time[idx]
  yevent <- survival_cohort$OS_event[idx]
  xb <- x[idx, , drop = FALSE]
  p <- cox_score_p(xb, ytime, yevent)
  names(p) <- colnames(xb)
  use <- names(p)[is.finite(p) & p < 0.05]
  selected <- character()
  if (length(use) > 1L) {
    inner <- integer(length(idx))
    for (event_value in 0:1) {
      j <- which(yevent == event_value)
      inner[j] <- sample(rep_len(1:10, length(j)))
    }
    cv <- cv.glmnet(xb[, use, drop = FALSE], Surv(ytime, yevent),
                    family = "cox", alpha = 1, standardize = TRUE,
                    foldid = inner, type.measure = "deviance")
    beta <- as.matrix(coef(cv, s = "lambda.min"))[, 1L]
    selected <- names(beta)[is.finite(beta) & beta != 0]
  }
  genes <- sym[intersect(selected, names(sym))]
  records[[b]] <- data.table(
    subsample = b, n = length(idx), events = sum(yevent),
    n_screened = length(use), n_selected = length(selected),
    gene = final8, selected = final8 %in% genes
  )
  if (b %% 10L == 0L) {
    cat("Completed 50% subsamples:", b, "/", n_subsamples, "\n")
    fwrite(rbindlist(records[seq_len(b)]),
           "results/half_sample_selection_detail.tsv", sep = "\t")
    flush.console()
  }
}

res <- rbindlist(records)
summary_dt <- res[, .(subsamples = uniqueN(subsample),
                      n_subsamples_selected = sum(selected),
                      selection_frequency = mean(selected)), by = gene]
cat("Writing final 50% subsample outputs\n")
fwrite(x = res, file = "results/half_sample_selection_detail.tsv", sep = "\t")
fwrite(x = summary_dt, file = "results/half_sample_selection_summary.tsv",
       sep = "\t")
pool_summary <- data.table(candidate_pool_size = length(candidate_ids))
fwrite(x = pool_summary, file = "results/half_sample_candidate_pool.tsv",
       sep = "\t")
capture.output(sessionInfo(), file = "results/half_sample_session.txt")
