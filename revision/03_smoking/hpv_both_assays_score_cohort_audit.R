suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
src <- c("Integration_Manuscript/Figure5_SharedNetworkArchitecture/Figure5_MasterData.tsv",
         "Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv")
paths <- c(file.path(root, src), "../01_baseline/results/sample_manifest.tsv")
m <- fread(paths[1], na.strings="NA")
a <- fread(paths[2], na.strings="NA")
b <- fread(paths[3], na.strings=c("NA", ""))[sample_type=="Primary Tumor"]
x <- merge(m, a, by.x="PatientID", by.y="patient_id")
x <- merge(x, b[,.(patient_barcode, primary_site, remove_base_and_mixed, p16, ISH)],
           by.x="PatientID", by.y="patient_barcode")
x[smoking_bin.x=="", smoking_bin.x:=NA_character_]
x[, smoke:=factor(smoking_bin.x, levels=c("Never", "Ever"))]
x[, sex:=factor(Sex)]; x[, stage:=factor(Stage)]
x[, site:=factor(fcase(grepl("Base of tongue",primary_site),"base_tongue",
                       grepl("pharynx",primary_site),"mixed",
                       grepl("tongue",primary_site),"other_tongue",
                       grepl("Floor of mouth",primary_site),"floor",
                       default="other_oral"))]
state <- function(z) fifelse(z %in% c("Negative", "Positive"), z, NA_character_)
x[, both_interpretable := !is.na(state(p16)) & !is.na(state(ISH))]
x[, hpv_status := factor(fcase(p16=="Positive" | ISH=="Positive", "Positive",
                               p16=="Negative" | ISH=="Negative", "Negative",
                               default=NA_character_), levels=c("Negative", "Positive"))]
scores <- c("Exposure2", "Proliferation2", "Stromal5")
clinical <- c("Age", "sex", "stage", "site", "hpv_status")
eligible_base <- complete.cases(x[, c("smoke", clinical), with=FALSE]) & x$both_interpretable
eligible_oral_base <- eligible_base & (x$remove_base_and_mixed %in% TRUE)
cat("Both-assay complete clinical cohort:", sum(eligible_base),
    "ever/never:", sum(x$smoke[eligible_base]=="Ever"), "/", sum(x$smoke[eligible_base]=="Never"),
    "oral subset:", sum(eligible_oral_base),
    "ever/never:", sum(x$smoke[eligible_oral_base]=="Ever"), "/", sum(x$smoke[eligible_oral_base]=="Never"), "\n")
for (y in scores) {
  eligible_y <- eligible_base & !is.na(x[[y]])
  oral_y <- eligible_y & (x$remove_base_and_mixed %in% TRUE)
  cat(y, ": N=", sum(eligible_y), " (", sum(x$smoke[eligible_y]=="Ever"), "/",
      sum(x$smoke[eligible_y]=="Never"), "); strict oral N=", sum(oral_y), " (",
      sum(x$smoke[oral_y]=="Ever"), "/", sum(x$smoke[oral_y]=="Never"), ")\n", sep="")
}
eligible <- eligible_base & complete.cases(x[, ..scores])
eligible_oral <- eligible & (x$remove_base_and_mixed %in% TRUE)
dir.create("hpv_audit", showWarnings=FALSE)
counts_for <- function(keep, label) data.table(
  cohort=label, n=sum(keep), ever=sum(x$smoke[keep]=="Ever",na.rm=TRUE),
  never=sum(x$smoke[keep]=="Never",na.rm=TRUE))
diagnostic <- rbindlist(list(
  counts_for(x$both_interpretable & complete.cases(x[,c("smoke",clinical),with=FALSE]),
             "Both assays interpretable; clinical/model covariates complete"),
  counts_for(eligible, "Both assays interpretable; all three scores and model covariates complete"),
  counts_for(eligible_oral, "Both assays interpretable; all three scores/covariates complete; oral-site restricted")
))
for (y in scores) {
  diagnostic <- rbind(diagnostic, counts_for(eligible_base & !is.na(x[[y]]),
                                               paste(y, "available; both assays and model covariates complete")))
}
fwrite(diagnostic, "hpv_audit/both_assays_score_cohort_diagnostic.tsv", sep="\t")
fwrite(data.table(input=c(paths, "../01_baseline/results/TCGA_HNSC_HPV_cBioPortal.tsv"),
  md5=unname(tools::md5sum(c(paths, "../01_baseline/results/TCGA_HNSC_HPV_cBioPortal.tsv")))),
  "hpv_audit/both_assays_score_cohort_audit_input_hashes.tsv", sep="\t")
capture.output(sessionInfo(), file="hpv_audit/both_assays_score_sensitivity_session.txt")
print(diagnostic)
