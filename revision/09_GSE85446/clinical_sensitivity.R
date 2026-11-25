library(survival)
d <- read.delim("cohort_audit.tsv")
con <- gzfile("GSE85446_series_matrix.txt.gz","rt")
header <- readLines(con,n=300)
close(con)
get_values <- function(s) gsub('"','',strsplit(s,'\t',fixed=TRUE)[[1]][-1],fixed=TRUE)
gsm <- get_values(header[startsWith(header,"!Sample_geo_accession\t")])
stopifnot(!anyDuplicated(gsm),!anyDuplicated(d$sample),setequal(gsm,d$sample))
get_meta <- function(key) {
  rows <- lapply(header[startsWith(header,"!Sample_characteristics_ch1\t")],get_values)
  row <- Filter(function(v) all(startsWith(v,paste0(key,": "))),rows)
  stopifnot(length(row)==1L)
  sub("^[^:]+: ","",row[[1]])[match(d$sample,gsm)]
}
st <- get_meta("pstagetnm")
sex <- get_meta("gender")
stopifnot(all(st %in% c("I","II","III","IV")),all(sex %in% c("f","m")))
d$late <- as.integer(st %in% c("III","IV"))
d$male <- as.integer(sex=="m")
d$age10 <- as.numeric(get_meta("age"))/10
stopifnot(nrow(d)==66,sum(d$event)==32,all(complete.cases(d)),all(d$age10>0),
  max(abs(as.numeric(scale(d$risk))-d$z))<1e-9)
forms <- list(clinical=Surv(months,event)~age10+male+late,
              combined=Surv(months,event)~age10+male+late+z)
models <- lapply(forms,coxph,data=d,ties="efron",x=TRUE)
fit <- models$combined
tab <- summary(fit)$coefficients
ci <- exp(confint(fit))
ph <- cox.zph(fit)$table
save_tsv <- function(x,file) write.table(x,file,sep="\t",quote=FALSE,row.names=FALSE)
adj <- data.frame(n=nrow(d),deaths=sum(d$event),late_stage=sum(d$late),
  HR_per_SD=exp(coef(fit)[["z"]]),low95=ci["z",1],high95=ci["z",2],P=tab["z",5],
  PH_score_P=ph["z","p"],PH_global_P=ph["GLOBAL","p"])
save_tsv(adj,"Final8_clinical_adjustment.tsv")
cstat <- function(model,data) {
  lp <- predict(model,newdata=data,type="lp")
  concordance(Surv(months,event)~lp,data=data,reverse=TRUE)$concordance
}
apparent <- vapply(models,cstat,0.0,data=d)

set.seed(85446)
B <- 500L
opt <- replicate(B,{
  boot <- d[sample.int(nrow(d),replace=TRUE),]
  tryCatch({
    fs <- lapply(forms,coxph,data=boot,ties="efron")
    if(any(!is.finite(unlist(lapply(fs,coef))))) stop("singular fit")
    vapply(fs,cstat,0.0,data=boot)-vapply(fs,cstat,0.0,data=d)
  },warning=function(w) c(NA_real_,NA_real_),error=function(e) c(NA_real_,NA_real_))
})
ok <- colSums(!is.finite(opt))==0
stopifnot(sum(ok)>=0.95*B)
corrected <- apparent-rowMeans(opt[,ok,drop=FALSE])
cmp <- data.frame(model=c(names(models),"combined_minus_clinical"),
  C_apparent=c(apparent,unname(diff(apparent))),
  C_optimism_corrected=c(corrected,unname(diff(corrected))),
  bootstrap_valid=sum(ok),bootstrap_requested=B)
save_tsv(cmp,"Final8_clinical_comparison.tsv")
print(adj)
print(cmp)
