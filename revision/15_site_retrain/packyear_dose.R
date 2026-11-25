suppressPackageStartupMessages({
  library(data.table)
})

root <- ".."
fig5 <- fread(paste0(Sys.getenv("OSCC_DATA"), "/Integration_Manuscript/Figure5_SharedNetworkArchitecture/Figure5_MasterData.tsv"))
clin <- fread(paste0(Sys.getenv("OSCC_DATA"), "/Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv"), na.strings = c("NA", ""))
manifest <- fread(file.path(root, "01_baseline/results/sample_manifest.tsv"),
                  na.strings = c("NA", ""))

scores <- c("Exposure2", "Proliferation2", "Stromal5")
stopifnot(all(c("PatientID", scores) %in% names(fig5)))
stopifnot(all(c("patient_id", "pack_years_smoked", "smoking_bin", "Age", "Sex", "Stage", "alcohol_bin") %in% names(clin)))
stopifnot(all(c("patient_barcode", "named_oral", "code") %in% names(manifest)))

site <- manifest[sample_type == "Primary Tumor" & remove_base_and_mixed == TRUE,
                 .(patient_id = patient_barcode, named_oral, code)]
dat <- merge(fig5[, c("PatientID", scores), with = FALSE],
             clin[, .(patient_id, pack_years_smoked, smoking_bin, Age, Sex,
                      Stage, alcohol_bin)],
             by.x = "PatientID", by.y = "patient_id", all = FALSE)
dat <- merge(dat, site, by.x = "PatientID", by.y = "patient_id", all = FALSE)
stopifnot(!anyDuplicated(dat$PatientID))
dat <- dat[named_oral == TRUE & smoking_bin == "Ever" &
             is.finite(pack_years_smoked) & pack_years_smoked > 0]
dat[, log2_pack_years := log2(pack_years_smoked)]
dat[, Sex := factor(Sex)]
dat[, Stage := factor(Stage)]
dat[, site_code := factor(code)]
dat[, alcohol_bin := factor(alcohol_bin)]

fit_models <- function(d, model_label, include_alcohol = FALSE) {
  rbindlist(lapply(scores, function(y) {
    covars <- c("log2_pack_years", "Age", "Sex", "Stage", "site_code",
                if (include_alcohol) "alcohol_bin")
    f <- reformulate(covars, response = y)
    m <- lm(f, data = d)
    cf <- summary(m)$coefficients
    est <- cf["log2_pack_years", ]
    ci <- confint(m)["log2_pack_years", ]
    data.table(model = model_label, outcome = y, n = nrow(d),
               estimate_per_doubling = est[["Estimate"]],
               conf_low = ci[[1L]], conf_high = ci[[2L]],
               p_value = est[["Pr(>|t|)"]], r_squared = summary(m)$r.squared)
  }))
}

base <- dat[complete.cases(dat[, c(scores, "log2_pack_years", "Age", "Sex", "Stage", "site_code"), with = FALSE])]
alcohol <- dat[complete.cases(dat[, c(scores, "log2_pack_years", "Age", "Sex", "Stage", "site_code", "alcohol_bin"), with = FALSE])]
setorder(base, PatientID)
setorder(alcohol, PatientID)
stopifnot(!anyDuplicated(base$PatientID), !anyDuplicated(alcohol$PatientID))

out <- rbind(fit_models(base, "primary_complete_case"),
             fit_models(alcohol, "alcohol_adjusted_complete_case", include_alcohol = TRUE))
out[, p_fdr := p.adjust(p_value, method = "BH"), by = model]
outdir <- file.path(root, "16_packyear_dose")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(out, file.path(outdir, "packyear_dose_models.tsv"), sep = "\t")
fwrite(data.table(model = c("primary_complete_case", "alcohol_adjusted_complete_case"),
                  n = c(nrow(base), nrow(alcohol)),
                  pack_year_median = c(median(base$pack_years_smoked), median(alcohol$pack_years_smoked)),
                  pack_year_q1 = c(quantile(base$pack_years_smoked, .25), quantile(alcohol$pack_years_smoked, .25)),
                  pack_year_q3 = c(quantile(base$pack_years_smoked, .75), quantile(alcohol$pack_years_smoked, .75))),
       file.path(outdir, "cohort_summary.tsv"), sep = "\t")
fwrite(base, file.path(outdir, "primary_analysis_cohort.tsv"), sep = "\t")
fwrite(alcohol, file.path(outdir, "alcohol_sensitivity_cohort.tsv"), sep = "\t")
cat("primary_n=", nrow(base), " alcohol_sensitivity_n=", nrow(alcohol), "\n", sep = "")
print(out)
