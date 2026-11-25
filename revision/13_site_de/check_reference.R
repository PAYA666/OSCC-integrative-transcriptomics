library(data.table)
root <- Sys.getenv("OSCC_DATA",Sys.getenv("OSCC_DATA"))
ref <- fread(file.path(root,"TCGA_OSCC_Nature8_DESeq2_ALL.tsv"),na.strings="NA")
archive <- readRDS(file.path(root,"TCGA_OSCC_Nature8_DEG_results.rds"))
strict <- function(x) with(x,is.finite(padj)&padj<.05&abs(log2FoldChange)>1)
stopifnot(setequal(ref$gene_id[strict(ref)],archive$gene_id),!anyDuplicated(ref$gene_id))
files <- c(paired_32="../01_baseline/results/paired_tumor_normal_DESeq2.tsv",
           remove_base_and_mixed="remove_base_and_mixed_DE.tsv", named_oral="named_oral_DE.tsv")
check <- function(i) {
  x <- fread(files[i],na.strings="NA")
  idx <- match(x$gene_id,ref$gene_id)
  common <- !is.na(idx)&is.finite(x$log2FoldChange)&is.finite(ref$log2FoldChange[idx])
  shared <- which(common & strict(x) & strict(ref)[idx])
  data.table(subset=names(files)[i],tumors=c(32,246,92)[i],normals=c(32,17,3)[i],
    common_genes=sum(common),rho=cor(x$log2FoldChange[common],ref$log2FoldChange[idx[common]],method="spearman"),
    strict_subset=sum(strict(x)),strict_reference=sum(strict(ref)),common_strict=length(shared),
    same_direction=sum(sign(x$log2FoldChange[shared])==sign(ref$log2FoldChange[idx[shared]])))
}
ans <- rbindlist(lapply(seq_along(files),check))
fwrite(ans,"verified_sensitivities.tsv",sep="\t")
print(ans)
