suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- Sys.getenv("OSCC_DATA")
expr <- readRDS(file.path(root, "1jieshu/04_External_Validation/GSE41613_frozen8_expression.rds"))
clinical <- readRDS(file.path(root, "GSE41613_Final8Gene_ExternalValidation.rds"))
stopifnot(nrow(expr) == 8L, ncol(expr) == 97L,
          setequal(colnames(expr), clinical$sample_id))
clinical <- clinical[match(colnames(expr), clinical$sample_id), ]
stopifnot(identical(clinical$sample_id, colnames(expr)),
          sum(clinical$OS_event) == 51L)

results <- rbindlist(lapply(rownames(expr), function(gene) {
  z <- as.numeric(scale(expr[gene, ]))
  fit_data <- data.table(OS_time_days = clinical$OS_time_days,
                         OS_event = clinical$OS_event, expression_z = z)
  fit <- coxph(Surv(OS_time_days, OS_event) ~ expression_z,
               data = fit_data, ties = "efron")
  ci <- exp(confint(fit)["expression_z", ])
  data.table(
    gene = gene,
    HR_per_external_SD = exp(coef(fit)["expression_z"]),
    lower95 = ci[1],
    upper95 = ci[2],
    p = summary(fit)$coefficients["expression_z", "Pr(>|z|)"],
    events = sum(fit_data$OS_event),
    n = nrow(fit_data)
  )
}))

fwrite(results, "GSE41613_individual_gene_HR.tsv", sep = "\t")
capture.output(sessionInfo(), file = "session.txt")
