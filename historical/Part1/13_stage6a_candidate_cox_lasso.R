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

cat("\nCandidate genes before survival screening:\n")
cat(nrow(candidate_precox), "\n")

print(
  table(
    candidate_precox$module_color,
    candidate_precox$DEG_status
  )
)

candidate_ids <- candidate_precox$gene_id

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

stopifnot(
  nrow(clinical_surv) == nrow(expr_surv)
)

cat("\nComplete OS patients:\n")
cat(nrow(clinical_surv), "\n")

cat("\nEvents:\n")
print(table(clinical_surv$OS_event))

dead_idx <- which(
  clinical_surv$OS_event == 1
)

alive_idx <- which(
  clinical_surv$OS_event == 0
)

train_dead <- sample(
  dead_idx,
  size = floor(length(dead_idx) * 0.70)
)

train_alive <- sample(
  alive_idx,
  size = floor(length(alive_idx) * 0.70)
)

train_idx <- sort(
  c(
    train_dead,
    train_alive
  )
)

valid_idx <- setdiff(
  seq_len(nrow(clinical_surv)),
  train_idx
)

clinical_train <- clinical_surv[
  train_idx,
  ,
  drop = FALSE
]

clinical_valid <- clinical_surv[
  valid_idx,
  ,
  drop = FALSE
]

expr_train <- expr_surv[
  train_idx,
  ,
  drop = FALSE
]

expr_valid <- expr_surv[
  valid_idx,
  ,
  drop = FALSE
]

cat("\nTraining cohort:\n")
cat(nrow(clinical_train), "\n")
print(table(clinical_train$OS_event))

cat("\nValidation cohort:\n")
cat(nrow(clinical_valid), "\n")
print(table(clinical_valid$OS_event))

split_table <- data.frame(
  patient_barcode = clinical_surv$patient_barcode,
  cohort = ifelse(
    seq_len(nrow(clinical_surv)) %in% train_idx,
    "Training",
    "Validation"
  ),
  OS_time = clinical_surv$OS_time,
  OS_event = clinical_surv$OS_event,
  stringsAsFactors = FALSE
)

write.table(
  split_table,
  "TCGA_OSCC_Nature8_OS_train_validation_split.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  split_table,
  "TCGA_OSCC_Nature8_OS_train_validation_split.rds"
)

uni_list <- vector(
  "list",
  length(candidate_ids)
)

names(uni_list) <- candidate_ids

for (i in seq_along(candidate_ids)) {

  gene <- candidate_ids[i]

  tmp <- data.frame(
    OS_time = clinical_train$OS_time,
    OS_event = clinical_train$OS_event,
    expr = expr_train[, gene]
  )

  fit <- try(
    survival::coxph(
      survival::Surv(
        OS_time,
        OS_event
      ) ~ expr,
      data = tmp
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {

    uni_list[[i]] <- data.frame(
      gene_id = gene,
      beta = NA_real_,
      HR = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      P = NA_real_
    )

  } else {

    s <- summary(fit)

    uni_list[[i]] <- data.frame(
      gene_id = gene,
      beta = s$coefficients[1, "coef"],
      HR = s$coefficients[1, "exp(coef)"],
      lower95 = s$conf.int[1, "lower .95"],
      upper95 = s$conf.int[1, "upper .95"],
      P = s$coefficients[1, "Pr(>|z|)"]
    )
  }
}

uni_cox_train <- dplyr::bind_rows(
  uni_list
)

uni_cox_train$FDR <- p.adjust(
  uni_cox_train$P,
  method = "BH"
)

anno_small <- candidate_precox %>%
  dplyr::select(
    gene_id,
    gene_name,
    module_color,
    MM,
    log2FoldChange,
    DEG_status
  )

uni_cox_train <- dplyr::left_join(
  uni_cox_train,
  anno_small,
  by = "gene_id"
)

uni_cox_train <- uni_cox_train %>%
  dplyr::arrange(P)

uni_selected <- uni_cox_train %>%
  dplyr::filter(
    !is.na(P),
    P < 0.05
  )

cat("\nTraining univariate Cox P < 0.05:\n")
cat(nrow(uni_selected), "\n")

cat("\nTraining univariate Cox FDR < 0.05:\n")
cat(
  sum(
    uni_cox_train$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

lasso_ids <- uni_selected$gene_id

if (length(lasso_ids) < 2) {

  stop(
    "Too few genes after training univariate Cox for LASSO."
  )
}

x_train <- as.matrix(
  expr_train[
    ,
    lasso_ids,
    drop = FALSE
  ]
)

y_train <- survival::Surv(
  clinical_train$OS_time,
  clinical_train$OS_event
)

set.seed(20260813)

cvfit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "cox",
  alpha = 1,
  nfolds = 10,
  standardize = TRUE,
  type.measure = "deviance"
)

cat("\nlambda.min:\n")
print(cvfit$lambda.min)

cat("\nlambda.1se:\n")
print(cvfit$lambda.1se)

coef_min <- coef(
  cvfit,
  s = "lambda.min"
)

idx_min <- which(
  as.numeric(coef_min) != 0
)

genes_min <- rownames(coef_min)[
  idx_min
]

beta_min <- as.numeric(coef_min)[
  idx_min
]

lasso_min <- data.frame(
  gene_id = genes_min,
  coefficient = beta_min,
  stringsAsFactors = FALSE
)

lasso_min <- dplyr::left_join(
  lasso_min,
  anno_small,
  by = "gene_id"
)

coef_1se <- coef(
  cvfit,
  s = "lambda.1se"
)

idx_1se <- which(
  as.numeric(coef_1se) != 0
)

genes_1se <- rownames(coef_1se)[
  idx_1se
]

beta_1se <- as.numeric(coef_1se)[
  idx_1se
]

lasso_1se <- data.frame(
  gene_id = genes_1se,
  coefficient = beta_1se,
  stringsAsFactors = FALSE
)

lasso_1se <- dplyr::left_join(
  lasso_1se,
  anno_small,
  by = "gene_id"
)

pdf(
  "Figure_Nature8_OSCC_LASSO_Cox_CV.pdf",
  width = 7,
  height = 6
)

plot(cvfit)

title(
  "Cross-validated LASSO-Cox",
  line = 2.5
)

dev.off()

png(
  "Figure_Nature8_OSCC_LASSO_Cox_CV.png",
  width = 1800,
  height = 1500,
  res = 250
)

plot(cvfit)

title(
  "Cross-validated LASSO-Cox",
  line = 2.5
)

dev.off()

lasso_path <- glmnet(
  x = x_train,
  y = y_train,
  family = "cox",
  alpha = 1,
  standardize = TRUE
)

pdf(
  "Figure_Nature8_OSCC_LASSO_Cox_CoefficientPath.pdf",
  width = 8,
  height = 6
)

plot(
  lasso_path,
  xvar = "lambda",
  label = FALSE
)

abline(
  v = log(cvfit$lambda.min),
  lty = 2
)

abline(
  v = log(cvfit$lambda.1se),
  lty = 3
)

dev.off()

write.table(
  candidate_precox,
  "TCGA_OSCC_Nature8_PreCox_candidate_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  uni_cox_train,
  "TCGA_OSCC_Nature8_Training_Univariate_Cox_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  uni_selected,
  "TCGA_OSCC_Nature8_Training_Univariate_Cox_P05.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  lasso_min,
  "TCGA_OSCC_Nature8_LASSO_lambdaMIN_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  lasso_1se,
  "TCGA_OSCC_Nature8_LASSO_lambda1SE_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  cvfit,
  "TCGA_OSCC_Nature8_LASSO_Cox_cvfit.rds"
)

saveRDS(
  lasso_min,
  "TCGA_OSCC_Nature8_LASSO_lambdaMIN_genes.rds"
)

saveRDS(
  lasso_1se,
  "TCGA_OSCC_Nature8_LASSO_lambda1SE_genes.rds"
)

saveRDS(
  clinical_train,
  "TCGA_OSCC_Nature8_training_clinical.rds"
)

saveRDS(
  clinical_valid,
  "TCGA_OSCC_Nature8_validation_clinical.rds"
)

saveRDS(
  expr_train,
  "TCGA_OSCC_Nature8_training_expression.rds"
)

saveRDS(
  expr_valid,
  "TCGA_OSCC_Nature8_validation_expression.rds"
)

cat("\n========================================\n")

cat(
  "Pre-Cox candidate genes:",
  nrow(candidate_precox),
  "\n"
)

cat(
  "Training patients:",
  nrow(clinical_train),
  "\n"
)

cat(
  "Validation patients:",
  nrow(clinical_valid),
  "\n"
)

cat(
  "Training univariate Cox P < 0.05:",
  nrow(uni_selected),
  "\n"
)

cat(
  "Genes at lambda.min:",
  nrow(lasso_min),
  "\n"
)

cat(
  "Genes at lambda.1se:",
  nrow(lasso_1se),
  "\n"
)

cat("========================================\n")
