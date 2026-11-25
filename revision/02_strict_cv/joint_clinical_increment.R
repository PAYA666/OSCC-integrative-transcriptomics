suppressPackageStartupMessages({library(data.table);library(survival);library(timeROC)})
revision_root <- "../"
data_root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
manifest <- fread(file.path(revision_root,"01_baseline/results/sample_manifest.tsv"), na.strings="NA")
folds <- fread(file.path(revision_root,"01_baseline/results/outer_folds.tsv"))
v <- manifest[sample_type=="Primary Tumor" & !is.na(OS_event)]
v[,stage_group:=factor(ifelse(Stage<=2,"I-II","III-IV"),levels=c("I-II","III-IV"))]
v[,age10:=as.numeric(Age)/10]
v[,Sex:=factor(Sex)]
clinical <- v[complete.cases(v[,.(OS_time,OS_event,Age,Sex,stage_group)])]
stopifnot(nrow(v)==297L,nrow(clinical)==275L,sum(clinical$OS_event)==140L,!anyDuplicated(clinical$patient_barcode))
raw <- readRDS(file.path(data_root,"TCGA_OSCC_Nature8_raw_count_matrix.rds"))
raw <- raw[raw$gene_type=="protein_coding",]
cnt <- as.matrix(raw[,manifest$sample_barcode,drop=FALSE]);rownames(cnt)<-raw$gene_id;storage.mode(cnt)<-"integer"
stopifnot(nrow(cnt)==19962L,ncol(cnt)==373L)
trans_vst <- function(q, dispersion) {a<-dispersion["extraPois"];b<-dispersion["asymptDisp"];log2((1+a+2*b*q+2*sqrt(b*q*(1+a+b*q)))/(4*b))}
risk_from_cox <- function(fit,newdata,times) {
  lp <- as.numeric(predict(fit,newdata=newdata,type="lp",reference="zero"))
  bh <- basehaz(fit,centered=FALSE)
  h0 <- vapply(times,function(tt){j<-findInterval(tt,bh$time);if(j<1)0 else bh$hazard[j]},numeric(1))
  rr <- 1-exp(-outer(exp(pmin(lp,700)),h0))
  rr
}
get_c <- function(d,marker) {
  stopifnot(marker %in% names(d))
  z <- copy(d)
  z$.marker <- z[[marker]]
  survival::concordance(Surv(OS_time,OS_event)~.marker,data=z,reverse=TRUE)$concordance
}

nfolds <- as.integer(Sys.getenv("N_FOLDS","100"));stopifnot(nfolds>=1L,nfolds<=100L)
files <- sprintf("checkpoints/r%02df%02d.rds",rep(1:10,each=10),rep(1:10,10))[seq_len(nfolds)]
allpred <- vector("list",length(files));foldstats <- vector("list",length(files))
times <- c(365.25,1095.75,1826.25)
for(i in seq_along(files)) {
  ck <- readRDS(files[i]);hold <- ck$predictions$patient_barcode
  tr <- v[!patient_barcode %in% hold];te <- v[match(hold,patient_barcode)]
  stopifnot(!anyNA(te$patient_barcode),!any(te$patient_barcode %in% tr$patient_barcode))
  ids <- names(ck$beta);beta <- ck$beta
  if(length(beta)==0L){stop(sprintf("No selected genes: %s",files[i]))}
  stopifnot(all(ids %in% ck$training_genes),all(ids %in% rownames(cnt)))
  tr_counts <- cnt[ids,tr$sample_barcode,drop=FALSE]
  tr_sf <- ck$size_factors[tr$sample_barcode]
  te_counts <- cnt[ids,te$sample_barcode,drop=FALSE]
  te_sf <- ck$held_size_factors[te$sample_barcode]
  stopifnot(!anyNA(tr_sf),!anyNA(te_sf),all(tr_sf>0),all(te_sf>0))
  tr_vst <- t(trans_vst(sweep(tr_counts,2,tr_sf,"/"),ck$dispersion))
  te_vst <- t(trans_vst(sweep(te_counts,2,te_sf,"/"),ck$dispersion))
  tr_lp <- as.numeric(tr_vst[,ids,drop=FALSE]%*%beta)-ck$center
  te_lp <- as.numeric(te_vst[,ids,drop=FALSE]%*%beta)-ck$center
  saved <- ck$predictions$lp[match(te$patient_barcode,ck$predictions$patient_barcode)]
  stopifnot(max(abs(te_lp-saved))<1e-7)
  tr_scores <- data.table(patient_barcode=tr$patient_barcode,score=tr_lp)
  te_scores <- data.table(patient_barcode=te$patient_barcode,score=te_lp)
  dtr <- merge(clinical,tr_scores,by="patient_barcode")
  dte <- merge(clinical,te_scores,by="patient_barcode")
  if(nrow(dtr)<200L || nrow(dte)<15L)stop(sprintf("Unexpected complete-case fold sizes in %s",files[i]))
  fitc <- coxph(Surv(OS_time,OS_event)~I(Age/10)+factor(Sex)+stage_group,data=dtr,ties="efron",x=TRUE)
  fitj <- coxph(Surv(OS_time,OS_event)~I(Age/10)+factor(Sex)+stage_group+score,data=dtr,ties="efron",x=TRUE)
  dte[,`:=`(genomic_lp=te_lp[match(patient_barcode,te$patient_barcode)],clinical_lp=as.numeric(predict(fitc,newdata=dte,type="lp")),joint_lp=as.numeric(predict(fitj,newdata=dte,type="lp")))]
  rc <- risk_from_cox(fitc,dte,times);rj <- risk_from_cox(fitj,dte,times)
  for(j in seq_along(times)){dte[,(paste0("clinical_risk_",c("1y","3y","5y")[j])):=rc[,j]];dte[,(paste0("joint_risk_",c("1y","3y","5y")[j])):=rj[,j]]}
  dte[,`:=`(repeat_id=ck$repeat_id,fold_id=ck$fold_id)]
  allpred[[i]] <- dte
  foldstats[[i]] <- data.table(repeat_id=ck$repeat_id,fold_id=ck$fold_id,n_train=nrow(dtr),n_test=nrow(dte),events_train=sum(dtr$OS_event),events_test=sum(dte$OS_event),selected_genes=length(beta),max_abs_reconstructed_score_error=max(abs(te_lp-saved)),joint_score_beta=unname(coef(fitj)["score"]))
  if(i%%10L==0L)cat("completed",i,"/",length(files),"\n")
}
oof <- rbindlist(allpred,fill=TRUE);stats <- rbindlist(foldstats,fill=TRUE)
if(nfolds==100L){stopifnot(nrow(oof)==2750L,all(oof[,.N,by=.(patient_barcode,repeat_id)]$N==1L),max(stats$max_abs_reconstructed_score_error)<1e-7)}

repeat_metrics <- oof[,.(genomic_C=get_c(.SD,"genomic_lp"),clinical_C=get_c(.SD,"clinical_lp"),joint_C=get_c(.SD,"joint_lp"),n=.N,events=sum(OS_event)),by=repeat_id]
repeat_metrics[,`:=`(joint_minus_clinical=joint_C-clinical_C,joint_minus_genomic=joint_C-genomic_C)]
patient <- oof[,.(OS_time=first(OS_time),OS_event=first(OS_event),genomic_lp=mean(genomic_lp),clinical_lp=mean(clinical_lp),joint_lp=mean(joint_lp),clinical_risk_1y=mean(clinical_risk_1y),clinical_risk_3y=mean(clinical_risk_3y),clinical_risk_5y=mean(clinical_risk_5y),joint_risk_1y=mean(joint_risk_1y),joint_risk_3y=mean(joint_risk_3y),joint_risk_5y=mean(joint_risk_5y),n_repeats=.N),by=patient_barcode]
if(nfolds==100L)stopifnot(nrow(patient)==275L,all(patient$n_repeats==10L))
setorder(patient,patient_barcode)
point <- c(genomic=get_c(patient,"genomic_lp"),clinical=get_c(patient,"clinical_lp"),joint=get_c(patient,"joint_lp"))
set.seed(20261120L);B<-2000L;boot<-matrix(NA_real_,nrow=B,ncol=3,dimnames=list(NULL,c("genomic","clinical","joint")))
for(b in seq_len(B)){ii<-sample.int(nrow(patient),nrow(patient),replace=TRUE);d<-patient[ii];for(m in colnames(boot))boot[b,m]<-tryCatch(get_c(d,paste0(m,"_lp")),error=function(e)NA_real_)}
boot<-boot[complete.cases(boot),,drop=FALSE]
metrics<-data.table(model=names(point),C_index=as.numeric(point),CI_lower=sapply(names(point),function(m)quantile(boot[,m],.025,names=FALSE)),CI_upper=sapply(names(point),function(m)quantile(boot[,m],.975,names=FALSE)),n=nrow(patient),events=sum(patient$OS_event),bootstrap_replicates=nrow(boot))
for(m in c("joint","genomic")) {dname<-paste0(m,"_minus_clinical");dv<-boot[,m]-boot[,"clinical"];metrics[model==m,`:=`(delta_vs_clinical=point[m]-point["clinical"],delta_CI_lower=quantile(dv,.025,names=FALSE),delta_CI_upper=quantile(dv,.975,names=FALSE))]}

auc <- rbindlist(lapply(c("genomic","clinical","joint"),function(m){z<-timeROC(T=patient$OS_time,delta=patient$OS_event,marker=patient[[paste0(m,"_lp")]],cause=1,times=times,iid=TRUE);ci<-confint(z)$CI_AUC/100;data.table(model=m,years=c(1,3,5),AUC=as.numeric(z$AUC),lower95=pmax(0,as.numeric(ci[,1])),upper95=pmin(1,as.numeric(ci[,2])))}))

calibrate <- function(model,horizon,col){
  d <- copy(patient)
  d <- d[is.finite(get(col)) & get(col)>=0 & get(col)<=1]
  if(nrow(d)<20L)stop("Too few patients for descriptive calibration")
  d[,risk:=get(col)]

  d[,risk_group:=pmin(5L,ceiling(data.table::frank(risk,ties.method="first")/.N*5L))]
  d[,{
    sf <- survfit(Surv(OS_time,OS_event)~1,data=.SD)
    obs <- summary(sf,times=horizon*365.25,extend=TRUE)$surv
    .(model=model,horizon_years=horizon,group=first(risk_group),n=.N,events=sum(OS_event),
      mean_predicted_risk=mean(risk),observed_KM_risk=1-as.numeric(obs),
      median_predicted_risk=median(risk))
  },by=risk_group]
}

calibration <- rbindlist(lapply(c("clinical","joint"),function(m) rbindlist(lapply(seq_along(times),function(j) calibrate(m,j,paste0(m,"_risk_",c("1y","3y","5y")[j]))))))
outdir <- "results"
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
fwrite(stats,file.path(outdir,sprintf("joint_clinical_fold_metrics_%03d_folds.tsv",nfolds)),sep="\t")
fwrite(repeat_metrics,file.path(outdir,sprintf("joint_clinical_repeat_metrics_%03d_folds.tsv",nfolds)),sep="\t")
fwrite(metrics,file.path(outdir,sprintf("joint_clinical_cindex_bootstrap_%03d_folds.tsv",nfolds)),sep="\t")
fwrite(auc,file.path(outdir,sprintf("joint_clinical_time_auc_%03d_folds.tsv",nfolds)),sep="\t")
fwrite(calibration,file.path(outdir,sprintf("joint_clinical_calibration_%03d_folds.tsv",nfolds)),sep="\t")
saveRDS(list(predictions=oof,patient_level=patient,fold_metrics=stats,repeat_metrics=repeat_metrics,cindex=metrics,time_auc=auc,calibration=calibration,bootstrap_replicates=B,bootstrap_scope="Patient bootstrap conditional on fixed repeated out-of-fold predictions; does not repeat model fitting."),file.path(outdir,sprintf("joint_clinical_increment_%03d_folds.rds",nfolds)))
cat("Completed",nfolds,"outer folds; patient-level N/events:",nrow(patient),sum(patient$OS_event),"\n")
print(metrics)
print(auc)
