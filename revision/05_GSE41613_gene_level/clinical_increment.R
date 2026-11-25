suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- Sys.getenv("OSCC_DATA")
pheno <- as.data.table(readRDS(file.path(root, "GSE41613_phenotype_local.rds")))
clinical <- as.data.table(readRDS(file.path(root,
  "GSE41613_Final8Gene_ExternalValidation.rds")))
stopifnot(nrow(pheno) == 97L, nrow(clinical) == 97L,
          sum(clinical$OS_event) == 51L)

get_characteristic <- function(x, prefix) {
  sub(paste0("^", prefix, ":\\s*"), "", x)
}
covariates <- pheno[, .(
  sample_id = geo_accession,
  age_group = get_characteristic(characteristics_ch1.2, "age"),
  sex = get_characteristic(characteristics_ch1.3, "Sex"),
  stage = get_characteristic(characteristics_ch1.1, "tumor stage"),
  treatment = get_characteristic(characteristics_ch1.6, "treatment")
)]
d <- merge(clinical, covariates, by = "sample_id", sort = FALSE)
stopifnot(nrow(d) == 97L, !anyDuplicated(d$sample_id))
d[, `:=`(age_group = factor(age_group), sex = factor(sex),
         stage_group = factor(stage), treatment = factor(treatment),
         risk_z = as.numeric(scale(risk_score)))]
stopifnot(all(is.finite(d$risk_z)),
          all(c("I/II", "III/IV") %in% levels(d$stage_group)))

clinical_formula <- Surv(OS_time_days, OS_event) ~ age_group + sex + stage_group
combined_formula <- update(clinical_formula, . ~ . + risk_z)
clinical_fit <- coxph(clinical_formula, data = d, ties = "efron", x = TRUE)
combined_fit <- coxph(combined_formula, data = d, ties = "efron", x = TRUE)

ci <- exp(confint(combined_fit)["risk_z", ])
lr <- 2 * (as.numeric(logLik(combined_fit)) - as.numeric(logLik(clinical_fit)))
clp <- predict(clinical_fit, type = "lp")
allp <- predict(combined_fit, type = "lp")
c_clinical <- concordance(Surv(OS_time_days, OS_event) ~ clp,
                          data = d, reverse = TRUE)$concordance
c_combined <- concordance(Surv(OS_time_days, OS_event) ~ allp,
                          data = d, reverse = TRUE)$concordance

set.seed(20261121L)
boot <- replicate(2000L, {
  idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
  b <- d[idx]
  tryCatch({
    m0 <- coxph(clinical_formula, data = b, ties = "efron")
    m1 <- coxph(combined_formula, data = b, ties = "efron")
    p0 <- predict(m0, newdata = b, type = "lp")
    p1 <- predict(m1, newdata = b, type = "lp")
    c0 <- concordance(Surv(OS_time_days, OS_event) ~ p0,
                      data = b, reverse = TRUE)$concordance
    c1 <- concordance(Surv(OS_time_days, OS_event) ~ p1,
                      data = b, reverse = TRUE)$concordance
    c(c0, c1, c1 - c0)
  }, error = function(e) c(NA_real_, NA_real_, NA_real_))
})
boot <- t(boot)
boot <- boot[complete.cases(boot), , drop = FALSE]

out <- data.table(
  model = c("Clinical only", "Clinical plus frozen Final8"),
  n = nrow(d), events = sum(d$OS_event),
  C_index_apparent = c(c_clinical, c_combined),
  C_index_bootstrap_median = c(median(boot[, 1]), median(boot[, 2])),
  C_index_bootstrap_lower95 = c(quantile(boot[, 1], .025),
                                quantile(boot[, 2], .025)),
  C_index_bootstrap_upper95 = c(quantile(boot[, 1], .975),
                                quantile(boot[, 2], .975)),
  HR_score_per_SD = c(NA_real_, exp(coef(combined_fit)["risk_z"])),
  HR_lower95 = c(NA_real_, ci[1]), HR_upper95 = c(NA_real_, ci[2]),
  score_p = c(NA_real_, summary(combined_fit)$coefficients["risk_z", "Pr(>|z|)"]),
  score_PH_p = c(NA_real_, cox.zph(combined_fit)$table["risk_z", "p"]),
  global_PH_p = c(cox.zph(clinical_fit)$table["GLOBAL", "p"],
                  cox.zph(combined_fit)$table["GLOBAL", "p"]),
  added_score_LR = c(NA_real_, lr),
  added_score_LR_p = c(NA_real_, pchisq(lr, df = 1L, lower.tail = FALSE)),
  bootstrap_replicates = nrow(boot),
  bootstrap_delta_C_median = c(NA_real_, median(boot[, 3])),
  bootstrap_delta_C_lower95 = c(NA_real_, quantile(boot[, 3], .025)),
  bootstrap_delta_C_upper95 = c(NA_real_, quantile(boot[, 3], .975))
)

fwrite(out, "GSE41613_clinical_increment.tsv", sep = "\t")
fwrite(d[, .(sample_id, OS_time_days, OS_event, risk_z, age_group, sex,
           stage_group, treatment)], "GSE41613_clinical_increment_data.tsv",
       sep = "\t")
capture.output(sessionInfo(), file = "GSE41613_clinical_increment_session.txt")
