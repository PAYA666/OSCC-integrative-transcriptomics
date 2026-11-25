source("score_final8.R")
root <- Sys.getenv("OSCC_DATA")
gettab <- function(name) read.delim(file.path(root,name))
raw <- gettab("TCGA_OSCC_Nature8_Final8Gene_CoxCoefficients.tsv")
std <- gettab("TCGA_OSCC_Final8Gene_STANDARDIZED_coefficients.tsv")
dat <- gettab("TCGA_OSCC_Nature8_Final8Gene_RiskScore.tsv")
v <- readRDS(file.path(root,"TCGA_OSCC_Nature8_VST_matrix.rds"))
tc <- colnames(v)[substr(colnames(v),14,15)=="01"]
patients <- substr(tc,1,12)
stopifnot(!anyDuplicated(patients), !anyDuplicated(dat$patient_barcode),
          nrow(dat)==297, all(dat$patient_barcode %in% patients),
          nrow(raw)==8, nrow(std)==8, !anyDuplicated(raw$gene_name),
          setequal(raw$gene_name,std$gene_name), all(raw$gene_id %in% rownames(v)))
raw <- raw[match(std$gene_name,raw$gene_name),]
x <- v[raw$gene_id,tc[match(dat$patient_barcode,patients)],drop=FALSE]
rownames(x) <- std$gene_name
ref <- data.frame(gene=std$gene_name,gene_id=raw$gene_id,center=rowMeans(x),
                  sd=apply(x,1,sd),beta=std$beta_z,row.names=NULL)
stopifnot(max(abs(raw$beta*ref$sd-ref$beta))<1e-10)
write.table(ref,"TCGA_reference.tsv",sep="\t",quote=FALSE,row.names=FALSE)
ref <- read.delim("TCGA_reference.tsv")
score <- function(y) score_final8(y,ref,"TCGA_Nature8_VST")
lp <- apply(x,2,score)
offset <- sum(raw$beta*ref$center)
err <- max(abs(lp+offset-dat$risk_score))
stopifnot(err<1e-10, abs(score(rev(x[,1]))-lp[1])<1e-12,
  inherits(try(score(x[-1,1]),silent=TRUE),"try-error"),
  inherits(try(score_final8(x[,1],ref,"microarray"),silent=TRUE),"try-error"))
report <- data.frame(n=length(lp),reference_raw_minus_standardized=offset,
                     max_abs_error_on_archived_scale=err,
                     standardized_score_median=median(lp))
write.table(report,"roundtrip_check.tsv",sep="\t",quote=FALSE,row.names=FALSE)
print(report)
