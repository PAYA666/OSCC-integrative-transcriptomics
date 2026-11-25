library(data.table)
library(survival)
library(timeROC)
f <- fread("../01_baseline/results/outer_folds.tsv")
d <- unique(fread("../01_baseline/results/sample_manifest.tsv",na.strings="NA")[
  sample_type=="Primary Tumor" & !is.na(OS_event),.(patient_barcode,OS_time,OS_event)])
ss <- list.files("checkpoints",pattern="^r[0-9]{2}f[0-9]{2}\\.rds$",full.names=TRUE)
stopifnot(length(ss)>0L)
x <- lapply(ss,readRDS)
stopifnot(all(vapply(x,function(a) identical(a$signature,x[[1]]$signature),logical(1))),
          all(vapply(x,function(a) !any(a$predictions$patient_barcode %in% a$discovery$patient_barcode),logical(1))))
p <- rbindlist(lapply(x,"[[","predictions"))
keys <- c("patient_barcode","repeat_id","fold_id")
stopifnot(!anyDuplicated(p[,..keys]),nrow(p[!f,on=keys])==0L,all(is.finite(p$lp)))
qc <- rbindlist(lapply(x,function(a) data.table(repeat_id=a$repeat_id,fold_id=a$fold_id,
    n_test=nrow(a$predictions),n_discovery=nrow(a$discovery),n_candidate=length(a$candidate),n_selected=length(a$beta))))
dir.create("results",showWarnings=FALSE)
fwrite(qc,"results/verified_folds.tsv",sep="\t")
cat("Verified folds:",length(x),"/100; scores:",nrow(p),"/2970\n")
if(length(x)!=100L) quit(save="no",status=0)
stopifnot(nrow(p)==2970L,nrow(d)==297L,all(p[,.N,by=patient_barcode]$N==10L))
z <- merge(d,p[,.(score=mean(lp)),by=patient_barcode],by="patient_barcode")
z[, score_z := as.numeric(scale(score))]
co <- concordance(Surv(OS_time,OS_event)~score_z,z,reverse=TRUE)
se <- sqrt(co$var)
fit <- coxph(Surv(OS_time,OS_event)~score_z,z)
score_ci <- exp(confint(fit)["score_z",])
auc <- timeROC(T=z$OS_time,delta=z$OS_event,marker=z$score_z,cause=1,
               times=c(365,1095,1825),iid=TRUE)
auc_ci <- confint(auc)$CI_AUC/100
repeat_c <- rbindlist(lapply(sort(unique(p$repeat_id)),function(r){
  q <- p[repeat_id==r,.(patient_barcode,score=lp)]
  stopifnot(nrow(q)==297L,!anyDuplicated(q$patient_barcode))
  zr <- merge(d,q,by="patient_barcode")
  stopifnot(nrow(zr)==297L)
  data.table(repeat_id=r,C=concordance(Surv(OS_time,OS_event)~score,zr,reverse=TRUE)$concordance)
}))
set.seed(20261120)
boot_c <- replicate(2000,{
  i <- sample.int(nrow(z),nrow(z),replace=TRUE)
  concordance(Surv(OS_time,OS_event)~score_z,z[i],reverse=TRUE)$concordance
})
fwrite(repeat_c,"results/repeat_cindex.tsv",sep="\t")
fwrite(p[,.(patient_barcode,repeat_id,score=lp)],
       "results/repeat_oof_predictions.tsv",sep="\t")
fwrite(data.table(years=c(1,3,5),AUC=unname(auc$AUC),
                  lower95=auc_ci[,1],upper95=auc_ci[,2]),
       "results/strict_pipeline_timeAUC.tsv",sep="\t")
fwrite(z,"results/patient_oof_scores.tsv",sep="\t")
fwrite(data.table(n=nrow(z),events=sum(z$OS_event),C=co$concordance,
  lower95_conditional=co$concordance-1.96*se,
  upper95_conditional=co$concordance+1.96*se,
  lower95_patient_bootstrap=unname(quantile(boot_c,.025)),
  upper95_patient_bootstrap=unname(quantile(boot_c,.975)),
  mean_repeat_C=mean(repeat_c$C),sd_repeat_C=sd(repeat_c$C),
  min_repeat_C=min(repeat_c$C),max_repeat_C=max(repeat_c$C),
  HR_per_SD=exp(coef(fit)["score_z"]),HR_lower95=score_ci[1],HR_upper95=score_ci[2],
  HR_p=summary(fit)$coefficients["score_z","Pr(>|z|)"],
  score_PH_p=cox.zph(fit)$table["score_z","p"],bootstrap_replicates=length(boot_c)),
  "results/strict_pipeline_cindex.tsv",sep="\t")
