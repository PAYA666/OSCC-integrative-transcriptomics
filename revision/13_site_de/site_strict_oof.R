library(data.table)
library(survival)
m <- fread("../01_baseline/results/sample_manifest.tsv", na.strings="NA")[sample_type=="Primary Tumor"]
p <- fread("../02_strict_cv/results/patient_oof_scores.tsv")
stopifnot(!anyDuplicated(m$patient_barcode), !anyDuplicated(p$patient_barcode))
d <- merge(p, m[,.(patient_barcode, OS_time, OS_event, remove_base_and_mixed, named_oral)],
           by=c("patient_barcode","OS_time","OS_event"))
stopifnot(nrow(d)==297L, sum(d$OS_event)==153L, sum(d$remove_base_and_mixed)==211L,
          sum(d$named_oral)==78L, !any(d$named_oral & !d$remove_base_and_mixed))
cstat <- function(z) concordance(Surv(OS_time, OS_event)~score_z, data=z, reverse=TRUE)$concordance
tab <- rbindlist(lapply(c("all","remove_base_and_mixed","named_oral"), function(k) {
    z <- if(k=="all") d else d[which(d[[k]])]
    set.seed(20261120)
    b <- replicate(2000, cstat(z[sample.int(nrow(z), replace=TRUE)]))
    data.table(subset=k, n=nrow(z), deaths=sum(z$OS_event), C=cstat(z),
               lower95=unname(quantile(b,.025)), upper95=unname(quantile(b,.975)))
}))
ref <- fread("../02_strict_cv/results/strict_pipeline_cindex.tsv")
stopifnot(abs(tab$C[1]-ref$C)<1e-10,all(is.finite(tab$lower95)))
fwrite(tab,"site_strict_oof.tsv",sep="\t")
print(tab)
