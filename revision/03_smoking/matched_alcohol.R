suppressPackageStartupMessages({library(DESeq2);library(data.table);library(fgsea)})
root <- Sys.getenv("OSCC_DATA")
dds0 <- readRDS(file.path(root,"Part2_TCGA_OSCC_Smoking_DESeq2_INPUT.rds"))
m <- as.data.table(as.data.frame(colData(dds0)))
man <- fread("../01_baseline/results/sample_manifest.tsv")
man <- man[sample_type=="Primary Tumor"]
stopifnot(ncol(dds0)==332L,nrow(dds0)==17389L,!anyDuplicated(m$patient_id),!anyDuplicated(man$patient_barcode))
idx <- match(m$patient_id,man$patient_barcode)
stopifnot(!anyNA(idx),identical(as.character(m$smoking_bin),as.character(man$smoking_bin[idx])))
m[,site:=factor(fcase(grepl("Base of tongue",man$primary_site[idx]),"base_tongue",
                     grepl("pharynx",man$primary_site[idx]),"mixed",
                     grepl("tongue",man$primary_site[idx]),"other_tongue",
                     grepl("Floor of mouth",man$primary_site[idx]),"floor",default="other_oral"))]
m[,sex:=factor(Sex)];m[,stage:=factor(Stage)];m[,age10:=as.numeric(as.character(Age))/10]
m[,smoke:=factor(smoking_bin,levels=c("Never","Ever"))]
ex <- fread(file.path(root,"Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv"),na.strings=c("NA",""))
stopifnot(!anyDuplicated(ex$patient_id),all(m$patient_id %in% ex$patient_id))
m[,alcohol:=factor(ex$alcohol_bin[match(patient_id,ex$patient_id)])]
vars <- c("age10","sex","stage","site","smoke","alcohol")
keep <- complete.cases(m[,..vars])
meta <- as.data.frame(droplevels(m[keep,..vars]))
ct <- counts(dds0)[,keep]
rownames(meta) <- colnames(ct)
stopifnot(identical(substr(colnames(ct),1,12),m$patient_id[keep]),!anyNA(meta))
dir.create("alcohol_matched",showWarnings=FALSE)
setwd("alcohol_matched")
f1 <- ~age10+sex+stage+site+smoke; f2 <- ~age10+sex+stage+site+alcohol+smoke
fwrite(as.data.table(meta)[,.N,by=.(smoke,alcohol)],"joint_counts.tsv",sep="\t")
X <- model.matrix(f2,meta)
stopifnot(qr(X)$rank==ncol(X),length(unique(meta$smoke))==2L)
main <- fread(file.path(root,"Part2_Smoking_DEG_Annotated.tsv"))
stopifnot(!anyDuplicated(main$gene_id),sum(main$padj<.05,na.rm=TRUE)==1275L)
gmt <- strsplit(readLines(file.path(root,"h.all.v2026.1.Hs.symbols.gmt")),"\t",fixed=TRUE)
pathways <- setNames(lapply(gmt,function(p) p[-(1:2)]),vapply(gmt,`[`,"",1L))
run <- function(f,label){
    message("Fitting ",label,"; N=",nrow(meta))
    d <- DESeqDataSetFromMatrix(ct,meta,f)
    fit <- DESeq(d,quiet=TRUE,parallel=FALSE)
    res <- results(fit, contrast=c("smoke","Ever","Never"),alpha=.05)
    rr <- as.data.table(as.data.frame(res),keep.rownames="gene_id")
    rr[,symbol:=main$gene_name[match(gene_id,main$gene_id)]]

    bad <- duplicated(rr$symbol)|duplicated(rr$symbol,fromLast=TRUE)
    rank <- rr[!bad & !is.na(symbol) & nzchar(symbol) & is.finite(stat)]
    setorder(rank,-stat,gene_id)
    ss <- setNames(rank$stat,rank$symbol)
    set.seed(20261121)
    g <- fgseaMultilevel(pathways=pathways,stats=ss,minSize=15,maxSize=500,nproc=1,eps=0)
    stopifnot(nrow(g)==50L)
    fwrite(rr,paste0("matched_genes_",label,".tsv"),sep="\t")
    g <- g[,.(pathway,size,NES,pval,padj)]
    fwrite(g,paste0("matched_hallmarks_",label,".tsv"),sep="\t")
    list(rr=rr,g=g)
}
u <- run(f1,"clinical");a <- run(f2,"clinical_plus_alcohol")
stopifnot(identical(u$rr$gene_id,a$rr$gene_id))
fold <- merge(u$rr,a$rr,by="gene_id",suffixes=c("_raw","_adj"))
joint <- fold[is.finite(log2FoldChange_raw)&is.finite(log2FoldChange_adj)]
common <- joint[padj_raw<.05 & padj_adj<.05]
summary <- data.table(n=nrow(meta),never=sum(meta$smoke=="Never"),ever=sum(meta$smoke=="Ever"),
    genes=nrow(ct),finite=nrow(joint),raw_fdr=sum(u$rr$padj<.05,na.rm=TRUE),adj_fdr=sum(a$rr$padj<.05,na.rm=TRUE),
    common=nrow(common),same_direction=sum(sign(common$log2FoldChange_raw)==sign(common$log2FoldChange_adj)),
    lfc_rho=cor(joint$log2FoldChange_raw,joint$log2FoldChange_adj,method="spearman"),
    raw_path_fdr=sum(u$g$padj<.05),adjusted_path_fdr=sum(a$g$padj<.05))
fwrite(summary,"matched_genome_summary.tsv",sep="\t")
fwrite(merge(u$g,a$g,by="pathway",suffixes=c("_raw","_adj")),"matched_hallmark_comparison.tsv",sep="\t")
capture.output(sessionInfo(),file="matched_genome_session.txt")
print(summary)
