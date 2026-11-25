suppressPackageStartupMessages({library(data.table);library(survival)})
root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
pheno <- as.data.table(readRDS(file.path(root, "GSE41613_phenotype_local.rds")))
clinical <- as.data.table(readRDS(file.path(root, "GSE41613_Final8Gene_ExternalValidation.rds")))
get_characteristic <- function(x, prefix) sub(paste0("^", prefix, ":\\s*"), "", x)
covars <- pheno[, .(
  sample_id=geo_accession,
  age_group=get_characteristic(characteristics_ch1.2, "age"),
  sex=get_characteristic(characteristics_ch1.3, "Sex"),
  stage_group=get_characteristic(characteristics_ch1.1, "tumor stage"),
  treatment=get_characteristic(characteristics_ch1.6, "treatment"),
  tissue_type=characteristics_ch1
)]
d <- merge(clinical, covars, by="sample_id", sort=FALSE)
stopifnot(nrow(d)==97L, !anyDuplicated(d$sample_id), sum(d$OS_event)==51L)
d[, risk_z:=as.numeric(scale(risk_score))]
d[, `:=`(age_group=factor(age_group), sex=factor(sex), stage_group=factor(stage_group),
         treatment=factor(treatment))]

d[tolower(as.character(treatment)) %in% c("unknown", "na", "not available"), treatment:=NA_character_]
d[as.character(treatment)=="multii-modality", treatment:="multi-modality"]

d <- d[!is.na(treatment) & complete.cases(risk_z, OS_time_days, OS_event, age_group, sex, stage_group)]
d[, treatment:=droplevels(treatment)]
cat("Treatment-complete cohort:", nrow(d), "events:", sum(d$OS_event),
    "treatment groups:", paste(names(table(d$treatment)),as.integer(table(d$treatment)),collapse="; "), "\n")
stopifnot(nrow(d)==96L, nlevels(droplevels(d$treatment))==2L)
f0 <- Surv(OS_time_days, OS_event) ~ risk_z + age_group + sex + stage_group
f1 <- update(f0, . ~ . + treatment)
m0 <- coxph(f0, data=d, ties="efron", x=TRUE)
m1 <- coxph(f1, data=d, ties="efron", x=TRUE)
extract <- function(model, label) {
  s <- summary(model); i <- match("risk_z", rownames(s$coefficients))
  ci <- exp(confint(model)["risk_z",])
  ph <- cox.zph(model, transform="km")
  data.table(model=label, n=nrow(d), events=sum(d$OS_event),
    treatment_counts=paste(names(table(d$treatment)), as.integer(table(d$treatment)), collapse="; "),
    HR_per_score_SD=exp(coef(model)["risk_z"]), lower95=ci[1], upper95=ci[2],
    p=s$coefficients[i,"Pr(>|z|)"], score_PH_p=ph$table["risk_z","p"],
    global_PH_p=ph$table["GLOBAL","p"])
}
out <- rbind(extract(m0,"Age/sex/stage adjusted; same treatment-complete cohort"),
             extract(m1,"Age/sex/stage/treatment adjusted"))
out[, treatment_source := "GSE41613 author-supplied phenotype; uni-modality vs multi-modality"]
fwrite(out,"GSE41613_treatment_adjusted_sensitivity.tsv",sep="\t")
inputs <- file.path(root,c("GSE41613_phenotype_local.rds","GSE41613_Final8Gene_ExternalValidation.rds"))
fwrite(data.table(input=inputs,md5=unname(tools::md5sum(inputs))),
       "GSE41613_treatment_adjusted_sensitivity_input_hashes.tsv",sep="\t")
capture.output(sessionInfo(),file="GSE41613_treatment_adjusted_sensitivity_session.txt")
print(out)
