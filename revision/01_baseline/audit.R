library(data.table)
library(survival)
root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
read <- function(name) fread(file.path(root, name), na.strings=c("NA",""))
dir.create("results", showWarnings=FALSE)
save_table <- function(x, name) fwrite(x, file.path("results",paste0(name,".tsv")),sep="\t",na="NA")

m <- read("TCGA_OSCC_Nature8_metadata.tsv")
c <- read("clinical.tsv")
v <- read("TCGA_OSCC_Nature8_Final8Gene_RiskScore.tsv")
s <- read("Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv")
hpv <- fread("results/TCGA_HNSC_HPV_cBioPortal.tsv",na.strings=c("NA",""))
dx <- unique(c[diagnoses.diagnosis_is_primary_disease == TRUE &
                 cases.submitter_id %in% m$patient_barcode,
               .(patient_barcode=cases.submitter_id,code=diagnoses.icd_10_code,
                 origin=diagnoses.tissue_or_organ_of_origin)])
stopifnot(!anyDuplicated(dx$patient_barcode), uniqueN(dx$patient_barcode)==341L,
          nrow(m)==373L, !anyDuplicated(m$sample_barcode),
          nrow(v)==297L,!anyDuplicated(v$patient_barcode),sum(v$OS_event)==153L,
          !anyDuplicated(s$patient_id))

named <- c("C02.1","C02.2","C03.0","C03.1","C03.9","C04.0","C04.9","C05.0","C06.0","C06.2")
dx[, named_oral := code %in% named]
dx[, remove_base_and_mixed := !code %in% c("C01","C14.8")]
x <- merge(m,dx,by="patient_barcode",all.x=TRUE)
x <- merge(x,v[,.(patient_barcode,OS_event,OS_time,risk_score)],by="patient_barcode",all.x=TRUE)
x <- merge(x,s[,.(patient_barcode=patient_id,smoking_bin,Age,Sex,Stage)],by="patient_barcode",all.x=TRUE)
x <- merge(x,hpv,by="patient_barcode",all.x=TRUE)
stopifnot(nrow(x)==373L, !anyNA(x$code), !anyDuplicated(hpv$patient_barcode))
save_table(x,"sample_manifest")
save_table(x[, .(samples=.N),by=.(code,origin,primary_site,named_oral,sample_type)], "site_codes")
summ <- function(d,label) {
 t <- d[sample_type=="Primary Tumor"]; n <- d[sample_type=="Solid Tissue Normal"]
 data.table(cohort=label,tumors=nrow(t),normals=nrow(n),paired=sum(n$patient_barcode %in% t$patient_barcode),
            survival=sum(!is.na(t$OS_event)),deaths=sum(t$OS_event,na.rm=TRUE),
            ever=sum(t$smoking_bin=="Ever",na.rm=TRUE),never=sum(t$smoking_bin=="Never",na.rm=TRUE))
}
save_table(rbind(summ(x,"submitted"),summ(x[remove_base_and_mixed==TRUE],"exclude_base_and_mixed"),
                 summ(x[named_oral==TRUE],"named_oral_subsites")), "cohort_counts")
save_table(x[sample_type=="Primary Tumor",
  .(n=.N),by=.(remove_base_and_mixed,named_oral,p16,ISH)],"hpv_by_site_definition")
save_table(x[sample_type=="Primary Tumor",
  .(p16_tested=sum(p16 %in% c("Positive","Negative")),
    p16_positive=sum(p16=="Positive",na.rm=TRUE),
    ISH_tested=sum(ISH %in% c("Positive","Negative")),
    ISH_positive=sum(ISH=="Positive",na.rm=TRUE)),by=.(remove_base_and_mixed,named_oral)],
  "hpv_coverage")

a <- read("TCGA_OSCC_Nature8_Stage6F_NestedCV_AllPredictions.tsv")
b <- read("TCGA_OSCC_Nature8_Stage6F_NestedCV_PatientPredictions.tsv")
stopifnot(nrow(a)==2970L,setequal(a$patient_barcode,v$patient_barcode),
          !anyDuplicated(a[,.(patient_barcode,repeat_id)]),
          identical(sort(unique(a$repeat_id)),1:10),
          all(a[,.(n=uniqueN(fold_id)),by=repeat_id]$n==10L))
j <- merge(a[,.(mean_prediction=mean(risk_score),sd_prediction=sd(risk_score)),by=patient_barcode],b,by="patient_barcode")
stopifnot(nrow(j)==297L,max(abs(j$mean_prediction-j$risk_score))<1e-10)
oof_c <- concordance(Surv(OS_time,OS_event)~mean_prediction,data=j,reverse=TRUE)$concordance
stopifnot(abs(oof_c-0.5695233)<1e-6)
save_table(data.table(n=nrow(j),events=sum(j$OS_event),predictions=nrow(a),
                      aggregation="patient arithmetic mean of 10 held-out scores",C=oof_c),"archived_oof_check")
save_table(a[,.(repeat_id,fold_id,patient_barcode)],"outer_folds")
qc <- a[,{
 ids <- patient_barcode
 train <- x[!patient_barcode %in% ids]
 stopifnot(!any(train$patient_barcode %in% ids))
 .(held_patients=.N,held_normals=sum(x$sample_type=="Solid Tissue Normal" & x$patient_barcode %in% ids),
   discovery_tumors=sum(train$sample_type=="Primary Tumor"),discovery_normals=sum(train$sample_type=="Solid Tissue Normal"))
},by=.(repeat_id,fold_id)]
save_table(qc,"outer_fold_sample_counts")

sub <- merge(x[sample_type=="Primary Tumor"],j[,.(patient_barcode,mean_prediction)],by="patient_barcode")
site_oof <- function(d,label){
 co <- concordance(Surv(OS_time,OS_event)~mean_prediction,data=d,reverse=TRUE)
 se <- sqrt(co$var)
 data.table(subset=label,n=nrow(d),events=sum(d$OS_event),C=co$concordance,
            lower95=co$concordance-1.96*se,upper95=co$concordance+1.96*se)
}
save_table(rbind(site_oof(sub,"submitted"),site_oof(sub[remove_base_and_mixed==TRUE],"exclude_base_and_mixed"),
                 site_oof(sub[named_oral==TRUE],"named_oral_subsites")),"site_archived_oof")
hpvneg <- x[sample_type=="Primary Tumor" & p16=="Negative" &
              !is.na(OS_event) & !is.na(risk_score)]
hpv_sensitivity <- function(z,label) {
 z[,risk_sd:=as.numeric(scale(risk_score))]
 f <- coxph(Surv(OS_time,OS_event)~risk_sd,data=z)
 ci <- exp(confint(f)["risk_sd",])
 data.table(subset=label,n=nrow(z),events=sum(z$OS_event),
   HR_per_SD=exp(coef(f)["risk_sd"]),lower95=ci[1],upper95=ci[2],
   p=summary(f)$coefficients["risk_sd","Pr(>|z|)"],
   C=concordance(f)$concordance)
}
save_table(rbind(hpv_sensitivity(hpvneg,"p16 negative"),
  hpv_sensitivity(hpvneg[remove_base_and_mixed==TRUE],"p16 negative, exclude base and mixed")),
  "p16_negative_Final8_sensitivity")

e <- read("GSE41613_Final8Gene_ExternalValidation.tsv")
p <- as.data.table(readRDS(file.path(root,"GSE41613_phenotype_local.rds")))
e <- merge(e,p[,.(sample_id=geo_accession,age_group=factor(`age:ch1`),sex=factor(`Sex:ch1`),
                   stage=factor(`tumor stage:ch1`),treatment=as.character(`treatment:ch1`))],by="sample_id")
stopifnot(!anyDuplicated(e$sample_id),nrow(e)==97L,sum(e$OS_event)==51L,all(e$OS_time_days>0))
e[, z:=as.numeric(scale(risk_score))]
stopifnot(max(abs(e$z-e$risk_z))<1e-8)
e[,treatment:=factor(fifelse(tolower(treatment)=="unknown",NA_character_,treatment))]
run <- function(d,label,adjust) {
 form <- reformulate(c("z",adjust),"Surv(OS_time_days,OS_event)")
 f <- coxph(form,d,x=TRUE,na.action=na.fail)
 ci <- exp(confint(f)["z",]); zz <- cox.zph(f)$table
 nullfit <- if(length(adjust)) coxph(reformulate(adjust,"Surv(OS_time_days,OS_event)"),data=d) else NULL
 lr <- if(length(adjust)) 2*as.numeric(logLik(f)-logLik(nullfit)) else NA_real_
 data.table(model=label,n=nrow(d),events=sum(d$OS_event),HR_per_external_SD=exp(coef(f)[["z"]]),
            lower95=ci[1],upper95=ci[2],score_PH_p=zz["z","p"],global_PH_p=zz["GLOBAL","p"],
            clinical_C_apparent=if(length(adjust)) concordance(nullfit)$concordance else NA_real_,
            C=concordance(f)$concordance,LR_added_score=lr,p_LR=pchisq(lr,1,lower.tail=FALSE))
}
save_table(rbind(run(e,"frozen score alone",NULL),run(e,"clinical + frozen score",c("age_group","sex","stage")),
                 run(e[!is.na(treatment)],"clinical + treatment + frozen score",c("age_group","sex","stage","treatment"))),
           "GSE41613_continuous_score")
files <- c("TCGA_OSCC_Nature8_metadata.tsv","clinical.tsv","TCGA_OSCC_Nature8_Final8Gene_RiskScore.tsv",
           "Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv","TCGA_OSCC_Nature8_Stage6F_NestedCV_AllPredictions.tsv",
           "TCGA_OSCC_Nature8_Stage6F_NestedCV_PatientPredictions.tsv","GSE41613_Final8Gene_ExternalValidation.tsv",
           "GSE41613_phenotype_local.rds")
file_paths <- c(file.path(root,files),"results/TCGA_HNSC_HPV_cBioPortal.tsv")
save_table(data.table(file=c(files,"cBioPortal TCGA HPV status"),md5=unname(tools::md5sum(file_paths))),"input_hashes")
capture.output(sessionInfo(),file="results/session.txt")
