suppressPackageStartupMessages({
  library(affy)
  library(AnnotationDbi)
  library(data.table)
  library(hgu133plus2.db)
  library(survival)
  library(timeROC)
})

root <- "."
out <- file.path(root, "results")
dir.create(out, showWarnings = FALSE)
manifest <- fread(file.path(root, "clinical_manifest.tsv"), na.strings = c("", "NA", "?"))
manifest[, `:=`(OS_time = as.numeric(futime), OS_event = as.integer(alldeath))]
manifest <- manifest[is.finite(OS_time) & OS_time > 0 & OS_event %in% 0:1]
stopifnot(nrow(manifest) == 71L, sum(manifest$OS_event) == 40L,
          all(manifest$`oc vs. op` == "Oral cavity"),
          all(manifest$`tu-vs-nl` == "cancer"))

files <- file.path(root, "raw", manifest$filename)
stopifnot(all(file.exists(files)),
          all(file.info(files)$size == as.numeric(manifest$archive_size)))
setorder(manifest, sample_id)
files <- file.path(root, "raw", manifest$filename)
arrays <- ReadAffy(filenames = files, verbose = FALSE)
sampleNames(arrays) <- manifest$sample_id
eset <- rma(arrays, verbose = FALSE)
probe_expr <- exprs(eset)
stopifnot(identical(colnames(probe_expr), manifest$sample_id))

coefficients <- fread(paste0(Sys.getenv("OSCC_DATA"), "/TCGA_OSCC_Final8Gene_STANDARDIZED_coefficients.tsv"))
genes <- coefficients$gene_name
beta <- setNames(coefficients$beta_z, genes)
symbols <- mapIds(hgu133plus2.db, keys = rownames(probe_expr),
                  column = "SYMBOL", keytype = "PROBEID", multiVals = "filter")
valid <- !is.na(symbols) & nzchar(symbols)
probe_index <- split(which(valid), symbols[valid])
aggregate_probes <- function(method) {
  reducer <- switch(method, median = stats::median, mean = base::mean)
  t(vapply(probe_index, function(i) {
    apply(probe_expr[i, , drop = FALSE], 2, reducer)
  }, numeric(ncol(probe_expr))))
}

gene_median <- aggregate_probes("median")
gene_mean <- aggregate_probes("mean")
stopifnot(all(genes %in% rownames(gene_median)),
          all(genes %in% rownames(gene_mean)))
score_expression <- function(expression) {
  z <- t(scale(t(expression[genes, , drop = FALSE])))
  raw_score <- as.numeric(crossprod(beta, z))
  as.numeric(scale(raw_score))
}

score_median <- score_expression(gene_median)
score_mean <- score_expression(gene_mean)
data <- copy(manifest)
data[, `:=`(risk_z = score_median, risk_z_mean_probe = score_mean,
            age_years = suppressWarnings(as.numeric(`age@dx`)),
            sex = factor(gender),
            stage_group = factor(ifelse(as.integer(ajcc4) <= 2, "I-II", "III-IV")))]
stopifnot(nrow(data) == 71L, all(is.finite(data$risk_z)))

fit_score <- function(d, score_name, label) {
  d <- copy(d)
  d[, score := get(score_name)]
  model <- coxph(Surv(OS_time, OS_event) ~ score, data = d, ties = "efron")
  ci <- exp(confint(model)["score", ])
  zph <- cox.zph(model)
  concordance_result <- concordance(Surv(OS_time, OS_event) ~ score, data = d,
                                    reverse = TRUE)
  data.table(
    analysis = label, n = nrow(d), events = sum(d$OS_event),
    HR_per_SD = exp(coef(model)["score"]), lower95_HR = ci[1],
    upper95_HR = ci[2], p = summary(model)$coefficients["score", "Pr(>|z|)"],
    C_index = concordance_result$concordance,
    C_lower95_conditional = concordance_result$concordance -
      1.96 * sqrt(concordance_result$var),
    C_upper95_conditional = concordance_result$concordance +
      1.96 * sqrt(concordance_result$var),
    PH_p = zph$table["score", "p"]
  )
}

primary <- fit_score(data, "risk_z", "Median probe summarization")
sensitivity_probe_mean <- fit_score(data, "risk_z_mean_probe", "Mean probe summarization")
adjusted <- data[complete.cases(data[, .(OS_time, OS_event, risk_z, age_years, sex, stage_group)])]
clinical_model <- coxph(Surv(OS_time, OS_event) ~ I(age_years / 10) + sex + stage_group,
                        data = adjusted, ties = "efron")
adjusted_model <- coxph(Surv(OS_time, OS_event) ~ risk_z + I(age_years / 10) + sex + stage_group,
                        data = adjusted, ties = "efron")
adjusted_ci <- exp(confint(adjusted_model)["risk_z", ])
lr_added <- 2 * (as.numeric(logLik(adjusted_model)) - as.numeric(logLik(clinical_model)))
clinical_lp <- predict(clinical_model, type = "lp")
combined_lp <- predict(adjusted_model, type = "lp")
adjusted_result <- data.table(
  analysis = "Score adjusted for age, sex, and AJCC stage",
  n = nrow(adjusted), events = sum(adjusted$OS_event),
  HR_per_SD = exp(coef(adjusted_model)["risk_z"]),
  lower95_HR = adjusted_ci[1], upper95_HR = adjusted_ci[2],
  p = summary(adjusted_model)$coefficients["risk_z", "Pr(>|z|)"],
  PH_p = cox.zph(adjusted_model)$table["risk_z", "p"],
  clinical_C_apparent = concordance(Surv(OS_time, OS_event) ~ clinical_lp,
                                    data = adjusted, reverse = TRUE)$concordance,
  clinical_plus_score_C_apparent = concordance(Surv(OS_time, OS_event) ~ combined_lp,
                                               data = adjusted, reverse = TRUE)$concordance,
  LR_added_score = lr_added, p_LR_added_score = pchisq(lr_added, 1, lower.tail = FALSE)
)

external_roc <- timeROC(T = data$OS_time, delta = data$OS_event,
                        marker = data$risk_z, cause = 1,
                        times = c(365, 1095, 1825), iid = TRUE)
external_auc_ci <- confint(external_roc)$CI_AUC / 100
external_auc <- data.table(years = c(1, 3, 5), AUC = unname(external_roc$AUC),
                           lower95 = external_auc_ci[, 1],
                           upper95 = external_auc_ci[, 2])

set.seed(20261120)
bootstrap_c <- replicate(2000, {
  i <- sample.int(nrow(data), nrow(data), replace = TRUE)
  d <- data[i]
  tryCatch(concordance(Surv(OS_time, OS_event) ~ risk_z, data = d,
                       reverse = TRUE)$concordance,
           error = function(e) NA_real_)
})
bootstrap_c <- bootstrap_c[is.finite(bootstrap_c)]
primary[, `:=`(C_lower95_bootstrap = unname(quantile(bootstrap_c, 0.025)),
               C_upper95_bootstrap = unname(quantile(bootstrap_c, 0.975)),
               bootstrap_replicates = length(bootstrap_c))]

probe_coverage <- data.table(gene = genes, n_probes = lengths(probe_index[genes]))
gene_output <- data.table(gene = genes,
                          as.data.frame(gene_median[genes, , drop = FALSE], check.names = FALSE))
score_output <- data[, .(sample_id, title, OS_time, OS_event, age_years, gender,
                         ajcc4, risk_z, risk_z_mean_probe)]
fwrite(probe_coverage, file.path(out, "Final8_probe_coverage.tsv"), sep = "\t")
fwrite(gene_output, file.path(out, "Final8_RMA_expression.tsv"), sep = "\t")
fwrite(score_output, file.path(out, "GSE42743_Final8_scores.tsv"), sep = "\t")
fwrite(rbindlist(list(primary, sensitivity_probe_mean), fill = TRUE),
       file.path(out, "GSE42743_external_validation.tsv"), sep = "\t")
fwrite(adjusted_result, file.path(out, "GSE42743_adjusted_validation.tsv"), sep = "\t")
fwrite(external_auc, file.path(out, "GSE42743_timeAUC.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "session.txt"))
