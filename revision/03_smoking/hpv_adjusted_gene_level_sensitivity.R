suppressPackageStartupMessages({library(DESeq2); library(data.table)})
root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
revision <- Sys.getenv("OSCC_REVISION", "..")
paths <- c(file.path(root,"Part2_TCGA_OSCC_Smoking_DESeq2_INPUT.rds"),
  file.path(root,"Integration_Manuscript/Figure5_SharedNetworkArchitecture/Figure5_MasterData.tsv"),
  file.path(root,"Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv"),
  file.path(revision,"01_baseline/results/sample_manifest.tsv"))
stopifnot(all(file.exists(paths)))
dds0 <- readRDS(paths[1]); m <- as.data.table(as.data.frame(colData(dds0)))
master <- fread(paths[2],na.strings="NA"); clin <- fread(paths[3],na.strings="NA")
man <- fread(paths[4],na.strings=c("NA",""))[sample_type=="Primary Tumor"]
stopifnot(!anyDuplicated(master$PatientID),!anyDuplicated(clin$patient_id),!anyDuplicated(man$patient_barcode),!anyDuplicated(m$patient_id))
x <- merge(master,clin,by.x="PatientID",by.y="patient_id")
x <- merge(x,man[,.(patient_barcode,primary_site,remove_base_and_mixed,p16,ISH)],by.x="PatientID",by.y="patient_barcode")
x[smoking_bin.x=="",smoking_bin.x:=NA_character_]
x[,smoke:=factor(smoking_bin.x,levels=c("Never","Ever"))]
x[,sex:=factor(Sex)]; x[,stage:=factor(Stage)]; x[,age10:=as.numeric(Age)/10]
x[,site:=factor(fcase(grepl("Base of tongue",primary_site),"base_tongue",
  grepl("pharynx",primary_site),"mixed",grepl("tongue",primary_site),"other_tongue",
  grepl("Floor of mouth",primary_site),"floor",default="other_oral"))]
state <- function(z) fifelse(z %in% c("Negative","Positive"),z,NA_character_)
x[,both_assays_disagree:=!is.na(state(p16)) & !is.na(state(ISH)) & state(p16)!=state(ISH)]
x[,hpv_status:=factor(fcase(p16=="Positive"|ISH=="Positive","Positive",
  p16=="Negative"|ISH=="Negative","Negative",default=NA_character_),levels=c("Negative","Positive"))]

i <- match(m$patient_id,x$PatientID); stopifnot(!anyNA(i))
z <- copy(x[i]); stopifnot(identical(as.character(m$smoking_bin),as.character(z$smoking_bin.x)))
base_vars <- c("smoke","age10","sex","stage","site")
full_vars <- c(base_vars,"hpv_status")
outdir <- "hpv_gene_level"
dir.create(outdir,showWarnings=FALSE)
run_set <- function(label, keep) {
  zz <- droplevels(z[keep]); n <- nrow(zz); zz[,age10:=age10-mean(age10)]
  stopifnot(n>0L, sum(zz$smoke=="Ever")>0L, sum(zz$smoke=="Never")>0L)
  ct <- counts(dds0)[,keep,drop=FALSE]
  rownames(zz) <- colnames(ct)
  zz <- as.data.frame(zz)
  rownames(zz) <- colnames(ct)

  designs <- list(clinical=as.formula("~ age10 + sex + stage + site + smoke"),
                  clinical_HPV=as.formula("~ age10 + sex + stage + site + hpv_status + smoke"))
  fits <- list(); summaries <- list()
  for (nm in names(designs)) {
    mm <- model.matrix(designs[[nm]],zz)
    if (qr(mm)$rank != ncol(mm) || nrow(mm)-ncol(mm)<5L) stop("Rank-deficient/overfit design: ",label," / ",nm)
    d <- DESeqDataSetFromMatrix(round(ct),zz,designs[[nm]])
    d <- DESeq(d,quiet=TRUE,parallel=FALSE)
    r <- results(d,contrast=c("smoke","Ever","Never"),alpha=.05)
    tab <- as.data.table(as.data.frame(r),keep.rownames="gene_id")
    tab[, model:=nm]
    setorder(tab,padj,pvalue)
    fwrite(tab,file.path(outdir,paste0(label,"_",nm,"_genes.tsv")),sep="\t",na="NA")
    fits[[nm]] <- tab
    summaries[[nm]] <- data.table(cohort=label,model=nm,n=n,ever=sum(zz$smoke=="Ever"),never=sum(zz$smoke=="Never"),
      hpv_negative=sum(zz$hpv_status=="Negative"),hpv_positive=sum(zz$hpv_status=="Positive"),
      n_genes=nrow(tab),finite=sum(is.finite(tab$log2FoldChange)&is.finite(tab$stat)),FDR05=sum(tab$padj<.05,na.rm=TRUE),
      residual_df=nrow(mm)-ncol(mm),design_rank=qr(mm)$rank,design_cols=ncol(mm))
  }
  a <- fits$clinical; b <- fits$clinical_HPV
  j <- merge(a[,.(gene_id,lfc_clin=log2FoldChange,padj_clin=padj)],b[,.(gene_id,lfc_hpv=log2FoldChange,padj_hpv=padj)],by="gene_id")
  ok <- is.finite(j$lfc_clin)&is.finite(j$lfc_hpv)
  f <- j[ok]
  cmp <- data.table(cohort=label,n_genes=nrow(f),spearman_lfc=cor(f$lfc_clin,f$lfc_hpv,method="spearman"),
    same_direction=sum(sign(f$lfc_clin)==sign(f$lfc_hpv)),clinical_FDR05=sum(f$padj_clin<.05,na.rm=TRUE),
    HPV_FDR05=sum(f$padj_hpv<.05,na.rm=TRUE),both_FDR05=sum(f$padj_clin<.05&f$padj_hpv<.05,na.rm=TRUE))
  fwrite(j,file.path(outdir,paste0(label,"_paired_effect_comparison.tsv")),sep="\t",na="NA")
  list(summary=rbindlist(summaries),comparison=cmp)
}
cc <- complete.cases(z[,..full_vars])
oral <- cc & (z$remove_base_and_mixed %in% TRUE)
cat("HPV-complete all-site cohort:",sum(cc),"Ever/Never:",sum(z$smoke[cc]=="Ever"),"/",sum(z$smoke[cc]=="Never"),"; strict oral:",sum(oral),"Ever/Never:",sum(z$smoke[oral]=="Ever"),"/",sum(z$smoke[oral]=="Never"),"\n")
r1 <- run_set("HPV_available_all_sites",cc)
r2 <- run_set("HPV_available_strict_oral",oral)
fwrite(rbindlist(list(r1$summary,r2$summary)),file.path(outdir,"model_cohort_summary.tsv"),sep="\t")
fwrite(rbindlist(list(r1$comparison,r2$comparison)),file.path(outdir,"effect_comparison_summary.tsv"),sep="\t")
fwrite(data.table(input=paths,md5=unname(tools::md5sum(paths))),file.path(outdir,"input_hashes.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
