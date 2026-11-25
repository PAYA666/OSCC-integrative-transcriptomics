suppressPackageStartupMessages({library(DESeq2); library(data.table)})
root <- Sys.getenv("OSCC_DATA",Sys.getenv("OSCC_DATA"))
m <- fread("../01_baseline/results/sample_manifest.tsv",na.strings="NA")
raw <- readRDS(file.path(root,"TCGA_OSCC_Nature8_raw_count_matrix.rds"))
raw <- raw[raw$gene_type=="protein_coding",]
ref <- fread(file.path(root,"TCGA_OSCC_Nature8_DESeq2_ALL.tsv"),na.strings="NA")
stopifnot(!anyDuplicated(ref$gene_id),!anyDuplicated(m$sample_barcode))
one <- function(flag) {
  s <- m[get(flag)==TRUE]
  group <- factor(ifelse(s$sample_type=="Primary Tumor","tumor","normal"),
                  levels=c("normal","tumor"))
  stopifnot(all(s$sample_type %in% c("Primary Tumor","Solid Tissue Normal")))
  cnt <- as.matrix(raw[,s$sample_barcode])
  rownames(cnt) <- raw$gene_id
  keep <- rowSums(cnt>=10)>=10
  cd <- data.frame(group,row.names=s$sample_barcode)
  dds <- DESeqDataSetFromMatrix(cnt[keep,,drop=FALSE],cd,~group)
  dds <- DESeq(dds,quiet=TRUE)
  r <- as.data.table(results(dds,contrast=c("group","tumor","normal"),alpha=.05),keep.rownames="gene_id")
  fwrite(r,paste0(flag,"_DE.tsv"),sep="\t")
  both <- merge(r[,.(gene_id,L=log2FoldChange,Q=padj)],
                ref[,.(gene_id,L0=log2FoldChange,Q0=padj)],by="gene_id")
  strict <- function(L,Q) is.finite(L)&is.finite(Q)&abs(L)>1&Q<.05
  both[,shared:=strict(L,Q)&strict(L0,Q0)]
  data.table(definition=flag,tumors=sum(group=="tumor"),normals=sum(group=="normal"),
    genes_tested=nrow(r),common_finite_fc=sum(is.finite(both$L)&is.finite(both$L0)),
    rho=cor(both$L,both$L0,use="complete.obs",method="spearman"),
    restricted_strict=sum(strict(r$log2FoldChange,r$padj)),shared_strict=sum(both$shared),
    shared_same_direction=both[shared==TRUE,sum(sign(L)==sign(L0))])
}
ans <- rbind(one("remove_base_and_mixed"),one("named_oral"))
stopifnot(identical(ans$tumors,c(246L,92L)),identical(ans$normals,c(17L,3L)))
fwrite(ans,"site_DE_summary.tsv",sep="\t")
capture.output(sessionInfo(),file="session.txt")
print(ans)
