suppressPackageStartupMessages({library(data.table); library(survival)})
meta <- fread("../01_baseline/results/sample_manifest.tsv",na.strings="NA")[remove_base_and_mixed==TRUE]
v <- meta[sample_type=="Primary Tumor" & !is.na(OS_event)]
folds <- fread("restricted_folds.tsv")
stopifnot(nrow(v)==211L, sum(v$OS_event)==103L, !anyDuplicated(v$patient_barcode))
allow_partial <- "--check-only" %in% commandArgs(TRUE)
paths <- sort(list.files("checkpoints",pattern="^r[0-9]{2}f[0-9]{2}\\.rds$",full.names=TRUE))
stopifnot(length(paths)>0L, allow_partial || length(paths)==100L)
stamp <- unname(tools::md5sum(c("run.R","../01_baseline/results/outer_folds.tsv",
  "../01_baseline/results/sample_manifest.tsv",file.path(Sys.getenv("OSCC_DATA",Sys.getenv("OSCC_DATA")),"TCGA_OSCC_Nature8_raw_count_matrix.rds"))))
p <- rbindlist(lapply(paths,function(path){
  a <- readRDS(path)
  held <- folds[repeat_id==a$repeat_id & fold_id==a$fold_id,patient_barcode]
  tr <- meta[!patient_barcode %in% held]
  z <- a$predictions
  stopifnot(identical(a$signature,stamp),
    basename(path)==sprintf("r%02df%02d.rds",a$repeat_id,a$fold_id),
    setequal(tr$sample_barcode,a$discovery$sample_barcode),
    !any(a$discovery$patient_barcode %in% held),
    setequal(setdiff(v$patient_barcode,held),a$inner$patient),
    !anyDuplicated(z$patient_barcode),setequal(held,z$patient_barcode),
    all(z$repeat_id==a$repeat_id & z$fold_id==a$fold_id),all(is.finite(z$lp)))
  z
}))
cat("Verified fits:",length(paths),"; held-out predictions:",nrow(p),"\n")
if(allow_partial) quit(status=0)
stopifnot(nrow(p)==2110L,nrow(unique(p[,.(repeat_id,patient_barcode)]))==2110L,
          all(p[,.N,by=patient_barcode]$N==10L))
z <- merge(v[,.(patient_barcode,OS_time,OS_event)],p[,.(score=mean(lp)),by=patient_barcode],by="patient_barcode")
stopifnot(nrow(z)==211L)
cc <- function(x) concordance(Surv(OS_time,OS_event)~score,data=x,reverse=TRUE)$concordance
set.seed(20261120)
b <- replicate(2000L,cc(z[sample.int(nrow(z),replace=TRUE)]))
ans <- data.table(n=nrow(z),deaths=sum(z$OS_event),repeats=10L,outer_fits=100L,
                  C=cc(z),lower95=unname(quantile(b,.025)),upper95=unname(quantile(b,.975)))
dir.create("results",showWarnings=FALSE)
fwrite(z,"results/patient_oof.tsv",sep="\t")
fwrite(ans,"results/restricted_pipeline_cindex.tsv",sep="\t")
print(ans)
