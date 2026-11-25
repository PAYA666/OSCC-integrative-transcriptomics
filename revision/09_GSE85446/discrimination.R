library(survival)
library(timeROC)
d <- read.delim("cohort_audit.tsv")
stopifnot(!anyDuplicated(d$sample), nrow(d)==66L, sum(d$event)==32L,
          all(complete.cases(d)), all(d$months>0))
set.seed(20261120)
cstat <- function(d) concordance(Surv(months,event) ~ z, data=d, reverse=TRUE)$concordance
b <- replicate(2000, cstat(d[sample.int(nrow(d), replace=TRUE), ]))
stopifnot(all(is.finite(b)))
h <- c(12,36,60)
r <- timeROC(T=d$months,delta=d$event,marker=d$z,cause=1,times=h,iid=TRUE)
ci <- confint(r)$CI_AUC / 100
stopifnot(all(is.finite(ci)),all(r$AUC >= ci[,1] & r$AUC <= ci[,2]))
cc <- data.frame(C=cstat(d),lower95=unname(quantile(b,.025)),upper95=unname(quantile(b,.975)),
                 valid=length(b))
au <- data.frame(months=h,AUC=unname(r$AUC),lower95=ci[,1],upper95=ci[,2],
                 deaths_by_horizon=sapply(h,function(t) sum(d$event==1 & d$months<=t)),
                 followed_beyond=sapply(h,function(t) sum(d$months>t)))
for (x in c("cc","au")) write.table(get(x),paste0("Final8_discrimination_",x,".tsv"),
                                    sep="\t",row.names=FALSE,quote=FALSE)
print(cc)
print(au)
