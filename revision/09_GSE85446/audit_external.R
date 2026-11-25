library(survival)
m <- readLines(gzfile("GSE85446_series_matrix.txt.gz"))
a <- readLines(gzfile("GPL6480.annot.gz"))
extract <- function(x, begin, end) {
  i <- which(x == begin); j <- which(x == end)
  stopifnot(length(i) == 1L, length(j) == 1L, j > i+1)
  read.delim(text = x[(i+1):(j-1)], check.names=FALSE, comment.char="",
             quote='"', na.strings=c("", "NA", "null"))
}
expr <- extract(m, "!series_matrix_table_begin", "!series_matrix_table_end")
ann <- extract(a, "!platform_table_begin", "!platform_table_end")
values <- function(line) gsub('"', '', strsplit(line, "\t",fixed=TRUE)[[1]][-1], fixed=TRUE)
field <- function(name) values(m[startsWith(m,paste0(name,"\t"))])
ch <- lapply(m[startsWith(m,"!Sample_characteristics_ch1\t")],values)
attribute <- function(key) {
  rows <- Filter(function(v) all(startsWith(v,paste0(key,": "))),ch)
  stopifnot(length(rows)==1)
  sub("^[^:]+: ","",rows[[1]])
}
id <- field("!Sample_geo_accession")
stopifnot(identical(id,colnames(expr)[-1]), !anyDuplicated(id),
          !anyDuplicated(field("!Sample_title")),!anyDuplicated(expr[[1]]),
          all(attribute("location")=="Oral Cavity"),
          all(attribute("tissue type")=="HPV-negative OSCC primary tumor"))
event <- as.integer(attribute("survival"))
times <- Filter(function(v) all(grepl("^(deathmonth|followuptime): ",v)),ch)
stopifnot(length(times)==1,all(event %in% 0:1), all(startsWith(times[[1]],"deathmonth: ")==(event==1)))
d <- data.frame(sample=id,event=event,months=as.numeric(sub("^[^:]+: ","",times[[1]])))
stopifnot(all(is.finite(d$months)),all(d$months>0))
b <- read.delim(paste0(Sys.getenv("OSCC_DATA"), "/TCGA_OSCC_Final8Gene_STANDARDIZED_coefficients.tsv"))

con <- gzfile("GSE85446_family.soft.gz","rt")
soft <- character()
repeat {
  chunk <- readLines(con,n=50000)
  if(!length(chunk)) break
  soft <- c(soft,grep("^(\\^SAMPLE = |!Sample_(source_name|label)_ch[12] = |#VALUE = )",chunk,value=TRUE))
}
close(con)
blocks <- split(soft,cumsum(startsWith(soft,"^SAMPLE = ")))
stopifnot(length(blocks)==length(id),
          setequal(vapply(blocks,function(v) sub("^\\^SAMPLE = ","",v[1]),""),id),
          all(vapply(blocks,function(v) all(c("!Sample_source_name_ch1 = HNSCC",
            "!Sample_label_ch1 = Cy3","!Sample_source_name_ch2 = reference",
            "!Sample_label_ch2 = Cy5","#VALUE = Normalized log2 ratios (Cy5/Cy3) representing reference/test") %in% v),logical(1))))
X <- -as.matrix(expr[,-1]); rownames(X) <- expr[[1]]
probes <- lapply(b$gene_name,function(g) intersect(ann$ID[which(ann[["Gene symbol"]]==g)],rownames(X)))
usable <- lapply(probes,function(p) p[vapply(p,function(i) all(is.finite(X[i,])) && sd(X[i,])>0, logical(1))])
tab <- data.frame(gene=b$gene_name,annotated_measured=lengths(probes),usable=lengths(usable),
                  probe_ids=vapply(usable,paste,collapse=";",FUN.VALUE=character(1)))
save_tsv <- function(x,file) write.table(x,file,sep="\t",row.names=FALSE,quote=FALSE)
save_tsv(tab,"Final8_probe_audit.tsv")
print(tab)
cat("Samples:",nrow(d),"events:",sum(d$event),"\n")
if(all(lengths(usable)>0)){
  genes <- sapply(usable,function(p) apply(X[p,,drop=FALSE],2,median))
  stopifnot(all(apply(genes,2,sd)>0))
  d$risk <- as.numeric(scale(genes) %*% b$beta_z)
  d$z <- as.numeric(scale(d$risk))
  fit <- coxph(Surv(months,event)~z,d)
  ci <- exp(confint(fit))
  cc <- concordance(Surv(months,event)~z,d,reverse=TRUE)
  se <- sqrt(cc$var)
  res <- data.frame(n=nrow(d),events=sum(event),HR_SD=unname(exp(coef(fit))),
    lower=ci[1],upper=ci[2],p=summary(fit)$coefficients[1,5],C=cc$concordance,
    C_lower=cc$concordance-qnorm(.975)*se,C_upper=cc$concordance+qnorm(.975)*se)
  save_tsv(res,"Final8_survival.tsv")
  print(res)
}
save_tsv(d,"cohort_audit.tsv")
