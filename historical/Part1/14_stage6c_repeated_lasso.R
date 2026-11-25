rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(glmnet)
library(dplyr)

set.seed(20260813)

gene_cox <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_candidateModules_geneCox_ALL.rds"
)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

clinical_wgcna <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_clinical_aligned.rds"
)

candidate_precox <- gene_cox %>%
  dplyr::filter(
    module_color %in% c("blue", "black"),
    abs(MM) >= 0.70,
    DEG_status %in% c("Up", "Down")
  ) %>%
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  )

candidate_ids <- candidate_precox$gene_id

cat("\nPre-Cox candidate genes:\n")
print(length(candidate_ids))

surv_keep <- (
  !is.na(clinical_wgcna$OS_time) &
  clinical_wgcna$OS_time > 0 &
  !is.na(clinical_wgcna$OS_event)
)

clinical_surv <- clinical_wgcna[
  surv_keep,
  ,
  drop = FALSE
]

expr_surv <- datExpr[
  surv_keep,
  candidate_ids,
  drop = FALSE
]

cat("\nPatients:\n")
print(nrow(clinical_surv))

cat("\nEvents:\n")
print(table(clinical_surv$OS_event))

x <- as.matrix(
  expr_surv
)

y <- survival::Surv(
  clinical_surv$OS_time,
  clinical_surv$OS_event
)

n_repeat <- 100
n_fold   <- 10

selection_matrix <- matrix(
  0,
  nrow = length(candidate_ids),
  ncol = n_repeat,
  dimnames = list(
    candidate_ids,
    paste0("rep", seq_len(n_repeat))
  )
)

coef_matrix <- matrix(
  0,
  nrow = length(candidate_ids),
  ncol = n_repeat,
  dimnames = list(
    candidate_ids,
    paste0("rep", seq_len(n_repeat))
  )
)

lambda_min_vec <- numeric(n_repeat)

for (r in seq_len(n_repeat)) {

  set.seed(
    20260813 + r
  )

  foldid <- rep(
    NA_integer_,
    nrow(clinical_surv)
  )

  dead_idx <- which(
    clinical_surv$OS_event == 1
  )

  alive_idx <- which(
    clinical_surv$OS_event == 0
  )

  foldid[dead_idx] <- sample(
    rep(
      seq_len(n_fold),
      length.out = length(dead_idx)
    )
  )

  foldid[alive_idx] <- sample(
    rep(
      seq_len(n_fold),
      length.out = length(alive_idx)
    )
  )

  cvfit <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    foldid = foldid,
    standardize = TRUE,
    type.measure = "deviance"
  )

  lambda_min_vec[r] <- cvfit$lambda.min

  beta <- as.numeric(
    coef(
      cvfit,
      s = "lambda.min"
    )
  )

  selection_matrix[, r] <- as.integer(
    beta != 0
  )

  coef_matrix[, r] <- beta

  if (r %% 10 == 0) {

    cat(
      "Completed",
      r,
      "/",
      n_repeat,
      "\n"
    )
  }
}

selection_count <- rowSums(
  selection_matrix
)

selection_frequency <- selection_count /
  n_repeat

mean_coef_all <- rowMeans(
  coef_matrix
)

mean_coef_selected <- sapply(
  seq_len(nrow(coef_matrix)),
  function(i) {

    z <- coef_matrix[i, ]

    z <- z[z != 0]

    if (length(z) == 0) {
      return(0)
    }

    mean(z)
  }
)

positive_frequency <- sapply(
  seq_len(nrow(coef_matrix)),
  function(i) {

    z <- coef_matrix[i, ]

    z <- z[z != 0]

    if (length(z) == 0) {
      return(NA_real_)
    }

    mean(z > 0)
  }
)

negative_frequency <- sapply(
  seq_len(nrow(coef_matrix)),
  function(i) {

    z <- coef_matrix[i, ]

    z <- z[z != 0]

    if (length(z) == 0) {
      return(NA_real_)
    }

    mean(z < 0)
  }
)

direction_consistency <- pmax(
  positive_frequency,
  negative_frequency,
  na.rm = TRUE
)

stability_table <- data.frame(
  gene_id = candidate_ids,
  selection_count = selection_count,
  selection_frequency = selection_frequency,
  mean_coef_all = mean_coef_all,
  mean_coef_selected = mean_coef_selected,
  positive_frequency = positive_frequency,
  negative_frequency = negative_frequency,
  direction_consistency = direction_consistency,
  stringsAsFactors = FALSE
)

anno <- candidate_precox %>%
  dplyr::select(
    gene_id,
    gene_name,
    module_color,
    MM,
    log2FoldChange,
    DEG_status
  )

stability_table <- dplyr::left_join(
  stability_table,
  anno,
  by = "gene_id"
)

stability_table <- stability_table %>%
  dplyr::arrange(
    dplyr::desc(selection_frequency),
    dplyr::desc(direction_consistency)
  )

write.table(
  stability_table,
  "TCGA_OSCC_Nature8_RepeatedLASSO_Stability.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  stability_table,
  "TCGA_OSCC_Nature8_RepeatedLASSO_Stability.rds"
)

saveRDS(
  selection_matrix,
  "TCGA_OSCC_Nature8_RepeatedLASSO_SelectionMatrix.rds"
)

saveRDS(
  coef_matrix,
  "TCGA_OSCC_Nature8_RepeatedLASSO_CoefficientMatrix.rds"
)

saveRDS(
  lambda_min_vec,
  "TCGA_OSCC_Nature8_RepeatedLASSO_lambdaMIN.rds"
)

pdf(
  "Figure_Nature8_OSCC_RepeatedLASSO_SelectionFrequency.pdf",
  width = 8,
  height = 6
)

hist(
  stability_table$selection_frequency,
  breaks = 20,
  xlab = "Selection frequency",
  main = "Repeated LASSO-Cox selection frequency"
)

dev.off()

top_plot <- stability_table %>%
  dplyr::slice_head(
    n = min(20, nrow(stability_table))
  )

pdf(
  "Figure_Nature8_OSCC_RepeatedLASSO_Top20.pdf",
  width = 9,
  height = 7
)

barplot(
  rev(top_plot$selection_frequency),
  names.arg = rev(top_plot$gene_name),
  horiz = TRUE,
  las = 1,
  xlab = "Selection frequency",
  main = "Top stable genes in repeated LASSO-Cox"
)

dev.off()

cat("\n========================================\n")

cat(
  "Candidates:",
  nrow(stability_table),
  "\n"
)

cat(
  "Selection frequency >= 0.8:",
  sum(
    stability_table$selection_frequency >= 0.8
  ),
  "\n"
)

cat(
  "Selection frequency >= 0.6:",
  sum(
    stability_table$selection_frequency >= 0.6
  ),
  "\n"
)

cat(
  "Selection frequency >= 0.5:",
  sum(
    stability_table$selection_frequency >= 0.5
  ),
  "\n"
)

cat(
  "Selection frequency >= 0.3:",
  sum(
    stability_table$selection_frequency >= 0.3
  ),
  "\n"
)

cat("========================================\n")
