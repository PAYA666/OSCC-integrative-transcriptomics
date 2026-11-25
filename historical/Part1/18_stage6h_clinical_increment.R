rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(timeROC)
library(dplyr)

dat <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_RiskScore.rds"
)

dat$sex <- factor(
  dat$sex,
  levels = c("female", "male")
)

dat$stage4 <- factor(
  dat$stage,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IVA",
    "Stage IVB"
  )
)

dat$risk_z <- as.numeric(
  scale(dat$risk_score)
)

model_dat <- dat %>%
  dplyr::filter(
    !is.na(OS_time),
    OS_time > 0,
    !is.na(OS_event),
    !is.na(age),
    !is.na(sex),
    !is.na(stage4),
    !is.na(risk_z)
  )

cat("\nComplete cases for clinical comparison:\n")
cat(nrow(model_dat), "\n")

cat("\nEvents:\n")
print(table(model_dat$OS_event))

cat("\nStage distribution:\n")
print(table(model_dat$stage4))

uni_risk <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_z,
  data = model_dat
)

uni_age <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ age,
  data = model_dat
)

uni_sex <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ sex,
  data = model_dat
)

uni_stage <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ stage4,
  data = model_dat
)

fit_clinical <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ age + sex + stage4,
  data = model_dat,
  x = TRUE,
  y = TRUE
)

fit_combined <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ age + sex + stage4 + risk_z,
  data = model_dat,
  x = TRUE,
  y = TRUE
)

fit_risk <- survival::coxph(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_z,
  data = model_dat,
  x = TRUE,
  y = TRUE
)

lrt <- anova(
  fit_clinical,
  fit_combined,
  test = "LRT"
)

get_C <- function(fit, data) {

  lp <- predict(
    fit,
    newdata = data,
    type = "lp"
  )

  survival::concordance(
    survival::Surv(
      data$OS_time,
      data$OS_event
    ) ~ lp,
    reverse = TRUE
  )$concordance
}

C_risk <- get_C(
  fit_risk,
  model_dat
)

C_clinical <- get_C(
  fit_clinical,
  model_dat
)

C_combined <- get_C(
  fit_combined,
  model_dat
)

times <- c(
  365,
  1095,
  1825
)

lp_risk <- predict(
  fit_risk,
  type = "lp"
)

lp_clinical <- predict(
  fit_clinical,
  type = "lp"
)

lp_combined <- predict(
  fit_combined,
  type = "lp"
)

roc_risk <- timeROC::timeROC(
  T = model_dat$OS_time,
  delta = model_dat$OS_event,
  marker = lp_risk,
  cause = 1,
  weighting = "marginal",
  times = times,
  iid = TRUE
)

roc_clinical <- timeROC::timeROC(
  T = model_dat$OS_time,
  delta = model_dat$OS_event,
  marker = lp_clinical,
  cause = 1,
  weighting = "marginal",
  times = times,
  iid = TRUE
)

roc_combined <- timeROC::timeROC(
  T = model_dat$OS_time,
  delta = model_dat$OS_event,
  marker = lp_combined,
  cause = 1,
  weighting = "marginal",
  times = times,
  iid = TRUE
)

extract_cox <- function(fit) {

  s <- summary(fit)

  data.frame(
    Variable = rownames(s$coefficients),
    beta = s$coefficients[, "coef"],
    HR = s$coefficients[, "exp(coef)"],
    lower95 = s$conf.int[, "lower .95"],
    upper95 = s$conf.int[, "upper .95"],
    P = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
}

combined_table <- extract_cox(
  fit_combined
)

clinical_table <- extract_cox(
  fit_clinical
)

model_compare <- data.frame(

  Model = c(
    "8-gene risk",
    "Clinical",
    "Clinical + 8-gene risk"
  ),

  C_index = c(
    C_risk,
    C_clinical,
    C_combined
  ),

  AUC_1y = c(
    roc_risk$AUC[1],
    roc_clinical$AUC[1],
    roc_combined$AUC[1]
  ),

  AUC_3y = c(
    roc_risk$AUC[2],
    roc_clinical$AUC[2],
    roc_combined$AUC[2]
  ),

  AUC_5y = c(
    roc_risk$AUC[3],
    roc_clinical$AUC[3],
    roc_combined$AUC[3]
  )
)

write.table(
  combined_table,
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalMultivariableCox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  model_compare,
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalModelComparison.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  fit_combined,
  "TCGA_OSCC_Nature8_Final8Gene_CombinedClinicalModel.rds"
)

saveRDS(
  model_dat,
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalCompleteCases.rds"
)

cat("\n========================================\n")
cat("STAGE 6H\n")
cat("========================================\n")

cat(
  "Complete cases:",
  nrow(model_dat),
  "\n"
)

cat("\nCombined model:\n")
print(combined_table)

cat("\nLikelihood-ratio test:\n")
print(lrt)

cat("\nModel performance:\n")
print(model_compare)

cat("========================================\n")
