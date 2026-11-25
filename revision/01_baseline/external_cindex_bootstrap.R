library(data.table)
library(survival)

d <- fread(paste0(Sys.getenv("OSCC_DATA"), "/GSE41613_Final8Gene_ExternalValidation.tsv"))
stopifnot(nrow(d) == 97L, sum(d$OS_event) == 51L,
          all(is.finite(d$risk_z)), all(d$OS_time_days > 0))
co <- concordance(Surv(OS_time_days, OS_event) ~ risk_z, data = d,
                  reverse = TRUE)
fit <- coxph(Surv(OS_time_days, OS_event) ~ risk_z, data = d)
ci <- exp(confint(fit)["risk_z", ])

set.seed(20261120)
boot_c <- replicate(5000, {
  i <- sample.int(nrow(d), nrow(d), replace = TRUE)
  concordance(Surv(OS_time_days, OS_event) ~ risk_z, data = d[i],
              reverse = TRUE)$concordance
})

result <- data.table(
  n = nrow(d), events = sum(d$OS_event), C_index = co$concordance,
  lower95_conditional = co$concordance - 1.96 * sqrt(co$var),
  upper95_conditional = co$concordance + 1.96 * sqrt(co$var),
  lower95_patient_bootstrap = unname(quantile(boot_c, .025)),
  upper95_patient_bootstrap = unname(quantile(boot_c, .975)),
  HR_per_SD = exp(coef(fit)["risk_z"]), HR_lower95 = ci[1], HR_upper95 = ci[2],
  HR_p = summary(fit)$coefficients["risk_z", "Pr(>|z|)"],
  bootstrap_replicates = length(boot_c)
)
fwrite(result, "results/GSE41613_cindex_bootstrap.tsv", sep = "\t")
