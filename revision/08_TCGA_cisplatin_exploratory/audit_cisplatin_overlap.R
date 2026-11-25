suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))

data_root <- Sys.getenv("OSCC_DATA")
manifest <- fread("../01_baseline/results/sample_manifest.tsv")
scores <- fread(file.path(
  data_root,
  "Part3_Cisplatin/Step4_TCGA_WGCNA_Projection/Part3_TCGA_ResistanceAssociated_Scores.tsv"
))
clinical <- fread(file.path(data_root, "clinical.tsv"))

cohort <- manifest[
  sample_type == "Primary Tumor" & !is.na(OS_event),
  .(patient_barcode, sample_barcode, primary_site, named_oral,
    remove_base_and_mixed, OS_time, OS_event, Age, Sex, Stage, p16, ISH)
]
ultra7 <- scores[, .(sample = Sample, Ultra7_Score)]
ultra7[, patient_barcode := substr(sample, 1L, 12L)]

cisplatin <- unique(
  clinical[
    project.project_id == "TCGA-HNSC" &
      grepl("cisplatin", treatments.therapeutic_agents, ignore.case = TRUE),
    .(patient_barcode = cases.submitter_id,
      agent = treatments.therapeutic_agents,
      treatment_type = treatments.treatment_type,
      treatment_start = treatments.days_to_treatment_start)
  ],
  by = "patient_barcode"
)

dat <- merge(cohort, ultra7, by = "patient_barcode")
dat <- merge(dat, cisplatin[, .(patient_barcode, cisplatin_record = TRUE,
                                agent, treatment_type, treatment_start)],
             by = "patient_barcode", all.x = TRUE)
dat[is.na(cisplatin_record), cisplatin_record := FALSE]

cat("Prognostic cohort:", nrow(cohort), "patients\n")
cat("Matched to Ultra7 scores:", nrow(dat), "patients\n")
cat("Explicit cisplatin records in matched cohort:",
    dat[cisplatin_record == TRUE, .N], "patients; deaths:",
    dat[cisplatin_record == TRUE, sum(OS_event)], "\n")
cat("Within named oral-cavity sites:",
    dat[cisplatin_record == TRUE & named_oral == TRUE, .N], "patients; deaths:",
    dat[cisplatin_record == TRUE & named_oral == TRUE, sum(OS_event)], "\n")
cat("Within oral-cavity, p16-negative cases:",
    dat[cisplatin_record == TRUE & named_oral == TRUE & p16 == "Negative", .N],
    "patients; deaths:",
    dat[cisplatin_record == TRUE & named_oral == TRUE & p16 == "Negative",
        sum(OS_event)], "\n")

print(dat[, .(n = .N, deaths = sum(OS_event)),
          by = .(cisplatin_record, named_oral, p16)][order(-n)])

fit_group <- function(d, label) {
  if (nrow(d) < 10L || sum(d$OS_event) < 5L) {
    return(data.table(group = label, n = nrow(d), deaths = sum(d$OS_event),
                      HR_per_SD = NA_real_, lower95 = NA_real_,
                      upper95 = NA_real_, p = NA_real_, note = "Insufficient data"))
  }
  d[, score_z := as.numeric(scale(Ultra7_Score))]
  fit <- coxph(Surv(OS_time, OS_event) ~ score_z, data = d)
  data.table(group = label, n = nrow(d), deaths = sum(d$OS_event),
             HR_per_SD = exp(coef(fit)[["score_z"]]),
             lower95 = exp(confint(fit)["score_z", 1L]),
             upper95 = exp(confint(fit)["score_z", 2L]),
             p = summary(fit)$coefficients["score_z", "Pr(>|z|)"],
             note = "Unadjusted exploratory association; not a treatment interaction")
}

cis <- dat[cisplatin_record == TRUE]
cis_oral <- cis[named_oral == TRUE]
cis_p16neg_oral <- cis[named_oral == TRUE & p16 == "Negative"]
association <- rbindlist(list(
  fit_group(cis, "All cohort cases with explicit cisplatin record"),
  fit_group(cis_oral, "Named oral-cavity cases with explicit cisplatin record"),
  fit_group(cis_p16neg_oral,
            "p16-negative named oral-cavity cases with explicit cisplatin record")
), fill = TRUE)
print(association)
fwrite(association, "TCGA_Ultra7_cisplatin_survival_exploratory.tsv", sep = "\t")

fwrite(dat, "TCGA_Ultra7_cisplatin_cohort_audit.tsv", sep = "\t")
capture.output(sessionInfo(), file = "audit_session.txt")
