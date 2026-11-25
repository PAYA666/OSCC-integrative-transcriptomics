suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(timeROC)
})

manifest <- fread("../01_baseline/results/sample_manifest.tsv", na.strings = "NA")
clinical <- unique(manifest[
  sample_type == "Primary Tumor" & !is.na(OS_event),
  .(patient_barcode, OS_time, OS_event, Age, Sex, Stage)
])
clinical[, stage_group := factor(ifelse(Stage <= 2, "I-II", "III-IV"),
                                 levels = c("I-II", "III-IV"))]
clinical <- clinical[complete.cases(clinical[, .(OS_time, OS_event, Age, Sex, stage_group)])]
stopifnot(nrow(clinical) == 275L, sum(clinical$OS_event) == 140L,
          !anyDuplicated(clinical$patient_barcode))

folds <- fread("../01_baseline/results/outer_folds.tsv")
score <- fread("results/patient_oof_scores.tsv")
score <- merge(clinical[, .(patient_barcode)],
               score[, .(patient_barcode, score)], by = "patient_barcode")
stopifnot(nrow(score) == nrow(clinical))
score_repeat <- fread("results/repeat_oof_predictions.tsv")

fit_clinical_fold <- function(rid, fid) {
  held <- folds[repeat_id == rid & fold_id == fid, patient_barcode]
  test <- clinical[patient_barcode %in% held]
  train <- clinical[!patient_barcode %in% held]
  model <- coxph(Surv(OS_time, OS_event) ~ I(Age / 10) + factor(Sex) + stage_group,
                 data = train, ties = "efron")
  lp <- predict(model, newdata = test, type = "lp")
  data.table(patient_barcode = test$patient_barcode, repeat_id = rid,
             fold_id = fid, clinical_lp = as.numeric(lp))
}

clinical_oof <- rbindlist(lapply(seq_len(10), function(r) {
  r_folds <- folds[repeat_id == r, sort(unique(fold_id))]
  r_predictions <- rbindlist(lapply(r_folds, function(k) fit_clinical_fold(r, k)))
  stopifnot(nrow(r_predictions) == nrow(clinical),
            !anyDuplicated(r_predictions$patient_barcode))
  r_predictions
}))

clinical_average <- clinical_oof[, .(clinical_lp = mean(clinical_lp)),
                                 by = "patient_barcode"]
comparison <- Reduce(function(x, y) merge(x, y, by = "patient_barcode"),
                     list(clinical, score, clinical_average))
stopifnot(nrow(comparison) == 275L)

c_index <- function(data, predictor) {
  formula <- reformulate(predictor, response = "Surv(OS_time, OS_event)")
  concordance(formula, data = data, reverse = TRUE)$concordance
}

repeat_comparison <- rbindlist(lapply(seq_len(10), function(r) {
  q <- merge(clinical, clinical_oof[repeat_id == r], by = "patient_barcode")
  q <- merge(q, score_repeat[repeat_id == r,
                             .(patient_barcode, genomic_lp = score)],
             by = "patient_barcode")
  stopifnot(nrow(q) == nrow(clinical))
  data.table(repeat_id = r, clinical_C = c_index(q, "clinical_lp"),
             Final8_C = c_index(q, "genomic_lp"))
}))
repeat_comparison[, delta_C := Final8_C - clinical_C]

set.seed(20261120)
boot <- replicate(2000, {
  i <- sample.int(nrow(comparison), nrow(comparison), replace = TRUE)
  d <- comparison[i]
  c(clinical = c_index(d, "clinical_lp"), Final8 = c_index(d, "score"))
})
boot_delta <- boot["Final8", ] - boot["clinical", ]

summary <- data.table(
  n = nrow(comparison), events = sum(comparison$OS_event),
  clinical_C = c_index(comparison, "clinical_lp"),
  clinical_C_lower95_patient_bootstrap =
    unname(quantile(boot["clinical", ], .025)),
  clinical_C_upper95_patient_bootstrap =
    unname(quantile(boot["clinical", ], .975)),
  Final8_C = c_index(comparison, "score"),
  Final8_C_lower95_patient_bootstrap =
    unname(quantile(boot["Final8", ], .025)),
  Final8_C_upper95_patient_bootstrap =
    unname(quantile(boot["Final8", ], .975)),
  delta_C_Final8_minus_clinical = c_index(comparison, "score") -
    c_index(comparison, "clinical_lp"),
  delta_lower95_patient_bootstrap = unname(quantile(boot_delta, .025)),
  delta_upper95_patient_bootstrap = unname(quantile(boot_delta, .975)),
  mean_repeat_clinical_C = mean(repeat_comparison$clinical_C),
  sd_repeat_clinical_C = sd(repeat_comparison$clinical_C),
  mean_repeat_Final8_C = mean(repeat_comparison$Final8_C),
  sd_repeat_Final8_C = sd(repeat_comparison$Final8_C)
)

auc_times <- c(365, 1095, 1825)
auc_results <- rbindlist(lapply(c("clinical_lp", "score"), function(marker) {
  roc <- timeROC(T = comparison$OS_time, delta = comparison$OS_event,
                 marker = comparison[[marker]], cause = 1, times = auc_times,
                 iid = TRUE)
  ci <- confint(roc)$CI_AUC / 100
  data.table(model = marker, years = c(1, 3, 5), AUC = unname(roc$AUC),
             lower95 = ci[, 1], upper95 = ci[, 2])
}))

fwrite(clinical_oof, "results/clinical_oof_predictions.tsv", sep = "\t")
fwrite(repeat_comparison, "results/repeat_clinical_comparison.tsv", sep = "\t")
fwrite(summary, "results/clinical_vs_Final8_nested_cv.tsv", sep = "\t")
fwrite(auc_results, "results/clinical_vs_Final8_timeAUC.tsv", sep = "\t")
