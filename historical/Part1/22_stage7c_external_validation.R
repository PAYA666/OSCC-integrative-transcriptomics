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

coef8_old <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_CoxCoefficients.rds"
)

frozen_genes <- c(
  "CTSG",
  "P4HA1",
  "NKAIN1",
  "LRRC38",
  "SPSB4",
  "SPINK6",
  "TNFRSF4",
  "DEFB1"
)

gene_map <- coef8_old[
  match(
    frozen_genes,
    coef8_old$gene_name
  ),
  c("gene_name", "gene_id")
]

stopifnot(
  all(gene_map$gene_name == frozen_genes)
)

tcga_keep <- (
  !is.na(clinical_wgcna$OS_time) &
  clinical_wgcna$OS_time > 0 &
  !is.na(clinical_wgcna$OS_event)
)

tcga_clin <- clinical_wgcna[
  tcga_keep,
  ,
  drop = FALSE
]

tcga_expr <- datExpr[
  tcga_keep,
  gene_map$gene_id,
  drop = FALSE
]

colnames(tcga_expr) <- gene_map$gene_name

tcga_z <- scale(
  tcga_expr,
  center = TRUE,
  scale = TRUE
)

tcga_z <- as.matrix(tcga_z)

tcga_model_data <- data.frame(
  OS_time = tcga_clin$OS_time,
  OS_event = tcga_clin$OS_event,
  tcga_z,
  check.names = FALSE
)

formula_z <- as.formula(
  paste(
    "Surv(OS_time, OS_event) ~",
    paste(
      frozen_genes,
      collapse = " + "
    )
  )
)

fit_z <- survival::coxph(
  formula_z,
  data = tcga_model_data,
  x = TRUE,
  y = TRUE
)

beta_z <- stats::coef(
  fit_z
)

cat("\n========================================\n")
cat("Frozen standardized TCGA coefficients\n")
cat("========================================\n")

print(beta_z)

risk_tcga_z <- as.numeric(
  tcga_z[
    ,
    frozen_genes,
    drop = FALSE
  ] %*%
    beta_z[
      frozen_genes
    ]
)

tcga_cutoff <- median(
  risk_tcga_z,
  na.rm = TRUE
)

cat("\nTCGA training cutoff:\n")
print(tcga_cutoff)

geo_expr8 <- readRDS(
  "GSE41613_frozen8_expression.rds"
)

geo_clinical <- readRDS(
  "GSE41613_OS_clinical.rds"
)

geo_expr <- t(
  geo_expr8[
    frozen_genes,
    ,
    drop = FALSE
  ]
)

geo_clinical <- geo_clinical[
  match(
    rownames(geo_expr),
    geo_clinical$sample_id
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(geo_expr),
    geo_clinical$sample_id
  )
)

geo_z <- scale(
  geo_expr,
  center = TRUE,
  scale = TRUE
)

geo_z <- as.matrix(
  geo_z
)

geo_risk <- as.numeric(
  geo_z[
    ,
    frozen_genes,
    drop = FALSE
  ] %*%
    beta_z[
      frozen_genes
    ]
)

geo_validation <- geo_clinical

geo_validation$risk_score <- geo_risk

geo_validation$risk_group <- ifelse(
  geo_validation$risk_score >= tcga_cutoff,
  "High",
  "Low"
)

geo_validation$risk_group <- factor(
  geo_validation$risk_group,
  levels = c(
    "Low",
    "High"
  )
)

cat("\n========================================\n")
cat("GSE41613 risk groups\n")
cat("========================================\n")

print(
  table(
    geo_validation$risk_group
  )
)

geo_validation$risk_z <- as.numeric(
  scale(
    geo_validation$risk_score
  )
)

cox_geo <- survival::coxph(
  survival::Surv(
    OS_time_days,
    OS_event
  ) ~ risk_z,
  data = geo_validation
)

c_geo <- survival::concordance(
  survival::Surv(
    OS_time_days,
    OS_event
  ) ~ risk_score,
  data = geo_validation,
  reverse = TRUE
)

fit_geo_km <- survival::survfit(
  survival::Surv(
    OS_time_days,
    OS_event
  ) ~ risk_group,
  data = geo_validation
)

km_geo <- survminer::ggsurvplot(
  fit_geo_km,
  data = geo_validation,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title = "GSE41613 external validation",
  legend.title = "Risk group",
  legend.labs = c(
    "Low",
    "High"
  )
)

pdf(
  "Figure_GSE41613_Final8Gene_External_KM.pdf",
  width = 7,
  height = 7
)

print(km_geo)

dev.off()

roc_geo <- timeROC::timeROC(
  T = geo_validation$OS_time_days,
  delta = geo_validation$OS_event,
  marker = geo_validation$risk_score,
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
  "Figure_GSE41613_Final8Gene_External_timeROC.pdf",
  width = 7,
  height = 6
)

plot(
  roc_geo,
  time = 365,
  lwd = 2,
  title = FALSE
)

plot(
  roc_geo,
  time = 1095,
  add = TRUE,
  lwd = 2
)

plot(
  roc_geo,
  time = 1825,
  add = TRUE,
  lwd = 2
)

legend(
  "bottomright",
  legend = c(

    paste0(
      "1-year AUC = ",
      round(
        roc_geo$AUC[1],
        3
      )
    ),

    paste0(
      "3-year AUC = ",
      round(
        roc_geo$AUC[2],
        3
      )
    ),

    paste0(
      "5-year AUC = ",
      round(
        roc_geo$AUC[3],
        3
      )
    )
  ),

  lwd = 2
)

title(
  "GSE41613 external validation"
)

dev.off()

coef_z_table <- data.frame(
  gene_name = names(beta_z),
  beta_z = as.numeric(beta_z),
  stringsAsFactors = FALSE
)

write.table(
  coef_z_table,
  "TCGA_OSCC_Final8Gene_STANDARDIZED_coefficients.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  geo_validation,
  "GSE41613_Final8Gene_ExternalValidation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  fit_z,
  "TCGA_OSCC_Final8Gene_STANDARDIZED_CoxModel.rds"
)

saveRDS(
  beta_z,
  "TCGA_OSCC_Final8Gene_STANDARDIZED_beta.rds"
)

saveRDS(
  geo_validation,
  "GSE41613_Final8Gene_ExternalValidation.rds"
)

saveRDS(
  roc_geo,
  "GSE41613_Final8Gene_External_timeROC.rds"
)

cat("\n========================================\n")
cat("GSE41613 EXTERNAL VALIDATION\n")
cat("========================================\n")

cat(
  "Patients:",
  nrow(geo_validation),
  "\n"
)

cat(
  "Events:",
  sum(geo_validation$OS_event),
  "\n"
)

cat(
  "C-index:",
  c_geo$concordance,
  "\n"
)

cat(
  "1-year AUC:",
  roc_geo$AUC[1],
  "\n"
)

cat(
  "3-year AUC:",
  roc_geo$AUC[2],
  "\n"
)

cat(
  "5-year AUC:",
  roc_geo$AUC[3],
  "\n"
)

cat("========================================\n")
