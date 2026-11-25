rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(survminer)
library(timeROC)
library(dplyr)

set.seed(20260813)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

clinical_wgcna <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_clinical_aligned.rds"
)

nested_selection <- readRDS(
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_GeneSelection.rds"
)

stable_table <- nested_selection %>%
  dplyr::filter(
    selection_frequency >= 0.80
  ) %>%
  dplyr::arrange(
    dplyr::desc(selection_frequency)
  )

stable_ids <- stable_table$gene_id

cat("\nFrozen signature:\n")

print(
  stable_table[
    ,
    c(
      "gene_name",
      "module_color",
      "selection_frequency"
    )
  ]
)

stopifnot(
  length(stable_ids) == 8
)

keep <- (
  !is.na(clinical_wgcna$OS_time) &
  clinical_wgcna$OS_time > 0 &
  !is.na(clinical_wgcna$OS_event)
)

clinical_final <- clinical_wgcna[
  keep,
  ,
  drop = FALSE
]

expr_final <- datExpr[
  keep,
  stable_ids,
  drop = FALSE
]

stopifnot(
  nrow(clinical_final) == 297
)

model_data <- data.frame(
  OS_time = clinical_final$OS_time,
  OS_event = clinical_final$OS_event,
  expr_final,
  check.names = FALSE
)

rhs <- paste(
  paste0("`", stable_ids, "`"),
  collapse = " + "
)

formula8 <- as.formula(
  paste(
    "Surv(OS_time, OS_event) ~",
    rhs
  )
)

fit8 <- survival::coxph(
  formula8,
  data = model_data,
  x = TRUE,
  y = TRUE
)

s8 <- summary(fit8)

coef8 <- data.frame(
  gene_id = rownames(s8$coefficients),
  beta = s8$coefficients[, "coef"],
  HR = s8$coefficients[, "exp(coef)"],
  lower95 = s8$conf.int[, "lower .95"],
  upper95 = s8$conf.int[, "upper .95"],
  P = s8$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

coef8 <- dplyr::left_join(
  coef8,
  stable_table[
    ,
    c(
      "gene_id",
      "gene_name",
      "module_color",
      "selection_frequency"
    )
  ],
  by = "gene_id"
)

coef8 <- coef8[
  match(
    stable_ids,
    coef8$gene_id
  ),
  ,
  drop = FALSE
]

beta <- stats::coef(fit8)

risk_score <- as.numeric(
  as.matrix(
    expr_final[
      ,
      names(beta),
      drop = FALSE
    ]
  ) %*%
    beta
)

clinical_final$risk_score <- risk_score

risk_cutoff <- median(
  clinical_final$risk_score,
  na.rm = TRUE
)

clinical_final$risk_group <- ifelse(
  clinical_final$risk_score >= risk_cutoff,
  "High",
  "Low"
)

clinical_final$risk_group <- factor(
  clinical_final$risk_group,
  levels = c(
    "Low",
    "High"
  )
)

risk_cox <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_score,
  data = clinical_final
)

c_final <- survival::concordance(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_score,
  data = clinical_final,
  reverse = TRUE
)

fit_km <- survival::survfit(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_group,
  data = clinical_final
)

km_final <- survminer::ggsurvplot(
  fit_km,
  data = clinical_final,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title = "Final 8-gene prognostic signature",
  legend.title = "Risk group",
  legend.labs = c(
    "Low",
    "High"
  )
)

pdf(
  "Figure_Nature8_OSCC_Final8Gene_KM.pdf",
  width = 7,
  height = 7
)

print(km_final)

dev.off()

roc_final <- timeROC::timeROC(
  T = clinical_final$OS_time,
  delta = clinical_final$OS_event,
  marker = clinical_final$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(
    365,
    1095,
    1825
  ),
  iid = TRUE
)

pdf(
  "Figure_Nature8_OSCC_Final8Gene_timeROC.pdf",
  width = 7,
  height = 6
)

plot(
  roc_final,
  time = 365,
  lwd = 2,
  title = FALSE
)

plot(
  roc_final,
  time = 1095,
  add = TRUE,
  lwd = 2
)

plot(
  roc_final,
  time = 1825,
  add = TRUE,
  lwd = 2
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "1-year AUC = ",
      round(roc_final$AUC[1], 3)
    ),
    paste0(
      "3-year AUC = ",
      round(roc_final$AUC[2], 3)
    ),
    paste0(
      "5-year AUC = ",
      round(roc_final$AUC[3], 3)
    )
  ),
  lwd = 2
)

title(
  "Final 8-gene prognostic signature"
)

dev.off()

B <- 500

optimism <- rep(
  NA_real_,
  B
)

set.seed(20260813)

for (b in seq_len(B)) {

  idx <- sample(
    seq_len(nrow(model_data)),
    replace = TRUE
  )

  boot_data <- model_data[
    idx,
    ,
    drop = FALSE
  ]

  boot_fit <- try(
    survival::coxph(
      formula8,
      data = boot_data
    ),
    silent = TRUE
  )

  if (inherits(
    boot_fit,
    "try-error"
  )) {
    next
  }

  lp_boot <- predict(
    boot_fit,
    newdata = boot_data,
    type = "lp"
  )

  C_boot <- survival::concordance(
    survival::Surv(
      boot_data$OS_time,
      boot_data$OS_event
    ) ~ lp_boot,
    reverse = TRUE
  )$concordance

  lp_test <- predict(
    boot_fit,
    newdata = model_data,
    type = "lp"
  )

  C_test <- survival::concordance(
    survival::Surv(
      model_data$OS_time,
      model_data$OS_event
    ) ~ lp_test,
    reverse = TRUE
  )$concordance

  optimism[b] <- C_boot - C_test

  if (b %% 100 == 0) {

    cat(
      "Bootstrap",
      b,
      "/",
      B,
      "\n"
    )
  }
}

mean_optimism <- mean(
  optimism,
  na.rm = TRUE
)

bootstrap_C <- c_final$concordance -
  mean_optimism

gene_cor <- cor(
  expr_final,
  method = "spearman",
  use = "pairwise.complete.obs"
)

diag(gene_cor) <- NA

max_abs_cor <- max(
  abs(gene_cor),
  na.rm = TRUE
)

write.table(
  coef8,
  "TCGA_OSCC_Nature8_Final8Gene_CoxCoefficients.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  clinical_final,
  "TCGA_OSCC_Nature8_Final8Gene_RiskScore.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_cor,
  "TCGA_OSCC_Nature8_Final8Gene_GeneCorrelation.tsv",
  sep = "\t",
  quote = FALSE
)

saveRDS(
  fit8,
  "TCGA_OSCC_Nature8_Final8Gene_CoxModel.rds"
)

saveRDS(
  coef8,
  "TCGA_OSCC_Nature8_Final8Gene_CoxCoefficients.rds"
)

saveRDS(
  clinical_final,
  "TCGA_OSCC_Nature8_Final8Gene_RiskScore.rds"
)

saveRDS(
  risk_cutoff,
  "TCGA_OSCC_Nature8_Final8Gene_RiskCutoff.rds"
)

cat("\n========================================\n")
cat("FINAL 8-GENE MODEL\n")
cat("========================================\n")

cat(
  "Patients:",
  nrow(clinical_final),
  "\n"
)

cat(
  "Events:",
  sum(clinical_final$OS_event),
  "\n"
)

cat(
  "Apparent C-index:",
  c_final$concordance,
  "\n"
)

cat(
  "Bootstrap corrected C-index:",
  bootstrap_C,
  "\n"
)

cat(
  "Mean optimism:",
  mean_optimism,
  "\n"
)

cat(
  "1-year AUC:",
  roc_final$AUC[1],
  "\n"
)

cat(
  "3-year AUC:",
  roc_final$AUC[2],
  "\n"
)

cat(
  "5-year AUC:",
  roc_final$AUC[3],
  "\n"
)

cat(
  "Maximum absolute gene-gene Spearman correlation:",
  max_abs_cor,
  "\n"
)

cat("========================================\n")
