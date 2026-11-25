suppressPackageStartupMessages({library(affy);library(data.table);library(survival)})
src <- "../05_GSE42743"
m <- fread(file.path(src,"clinical_manifest.tsv"))
d <- fread(file.path(src,"results/GSE42743_Final8_scores.tsv"))
idx <- match(d$sample_id,m$sample_id)
stopifnot(nrow(d)==71,!anyDuplicated(d$sample_id),!anyNA(idx),
          all(d$OS_time==as.numeric(m$futime[idx])),all(d$OS_event==as.integer(m$alldeath[idx])))
files <- file.path(src,"raw",m$filename[idx])
stopifnot(all(file.info(files)$size==m$archive_size[idx]))
probes <- c("1552283_s_at","1554008_at","1568765_at","201108_s_at","201497_x_at",
            "205200_at","207517_at","209900_s_at","210797_s_at","222344_at",
            "230104_s_at","240000_at","242417_at")
beta <- c(.000443,.231391,.383316,.205039,-.206855,-.131867,-.008777,.171380,
          .216877,.199328,-.124634,.864964,-.438477)
arr <- ReadAffy(filenames=files,verbose=FALSE);sampleNames(arr)<-d$sample_id
ex <- exprs(rma(arr,verbose=FALSE))
stopifnot(identical(colnames(ex),d$sample_id),all(probes %in% rownames(ex)))
d[,old13:=as.numeric(scale(drop(crossprod(beta,ex[probes,,drop=FALSE]))))]
stopifnot(all(is.finite(d$old13)),sum(d$OS_event)==40)
cc <- function(x){
    old <- concordance(Surv(OS_time,OS_event)~old13,x,reverse=TRUE)$concordance
    new <- concordance(Surv(OS_time,OS_event)~risk_z,x,reverse=TRUE)$concordance
    c(published13=old,Final8=new,Final8_minus_published13=new-old)
}
obs <- cc(d)

set.seed(20261120);b <- replicate(2000,cc(d[sample.int(nrow(d),nrow(d),replace=TRUE)]))
stopifnot(all(is.finite(b)),abs(obs["Final8"]-.593)<.001)
tab <- data.table(estimand=names(obs),estimate=as.numeric(obs),
                  lower95=apply(b,1,quantile,.025),upper95=apply(b,1,quantile,.975),
                  n=nrow(d),deaths=sum(d$OS_event),bootstrap=2000)
fwrite(tab,"GSE42743_paired_signature_comparison.tsv",sep="\t")
fwrite(data.table(probe=probes,beta=beta),"published_table2_parameters.tsv",sep="\t")
fwrite(d[,.(sample_id,OS_time,OS_event,Final8_z=risk_z,published13_z=old13)],
       "GSE42743_score_audit.tsv",sep="\t")
capture.output(sessionInfo(),file="session.txt")
print(tab)
