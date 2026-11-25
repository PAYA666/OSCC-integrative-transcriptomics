rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(WGCNA)

options(stringsAsFactors = FALSE)

allowWGCNAThreads()

vst_mat <- readRDS(
  "TCGA_OSCC_Nature8_VST_matrix.rds"
)

meta <- readRDS(
  "TCGA_OSCC_Nature8_metadata.rds"
)

tumor_meta <- meta[
  meta$sample_type == "Primary Tumor",
  ,
  drop = FALSE
]

tumor_ids <- tumor_meta$sample_barcode

stopifnot(
  all(tumor_ids %in% colnames(vst_mat))
)

tumor_vst <- vst_mat[
  ,
  tumor_ids,
  drop = FALSE
]

cat("Tumor VST matrix:\n")
print(dim(tumor_vst))

gene_sd <- apply(
  tumor_vst,
  1,
  sd,
  na.rm = TRUE
)

tumor_vst <- tumor_vst[
  is.finite(gene_sd) &
    gene_sd > 0,
  ,
  drop = FALSE
]

cat("\nGenes after zero-variance removal:\n")
print(nrow(tumor_vst))

gene_mad <- apply(
  tumor_vst,
  1,
  mad,
  na.rm = TRUE
)

gene_mad <- sort(
  gene_mad,
  decreasing = TRUE
)

top_n <- min(
  5000,
  length(gene_mad)
)

top_genes <- names(
  gene_mad
)[1:top_n]

wgcna_expr <- tumor_vst[
  top_genes,
  ,
  drop = FALSE
]

cat("\nTop MAD genes selected:\n")
print(nrow(wgcna_expr))

datExpr <- t(wgcna_expr)

cat("\nWGCNA datExpr dimensions:\n")
print(dim(datExpr))

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

cat("\nAll samples/genes good:\n")
print(gsg$allOK)

if (!gsg$allOK) {

  if (sum(!gsg$goodGenes) > 0) {

    cat("\nBad genes:\n")

    print(
      colnames(datExpr)[
        !gsg$goodGenes
      ]
    )

  }

  if (sum(!gsg$goodSamples) > 0) {

    cat("\nBad samples:\n")

    print(
      rownames(datExpr)[
        !gsg$goodSamples
      ]
    )

  }

}

sampleTree <- hclust(
  dist(datExpr),
  method = "average"
)

pdf(
  "Figure_Nature8_OSCC_WGCNA_SampleClustering.pdf",
  width = 14,
  height = 7
)

plot(
  sampleTree,
  main = "Sample clustering for WGCNA",
  sub = "",
  xlab = "",
  cex.lab = 1.2,
  cex.axis = 1.2,
  cex.main = 1.4
)

dev.off()

png(
  "Figure_Nature8_OSCC_WGCNA_SampleClustering.png",
  width = 2400,
  height = 1200,
  res = 180
)

plot(
  sampleTree,
  main = "Sample clustering for WGCNA",
  sub = "",
  xlab = "",
  cex.lab = 1.2,
  cex.axis = 1.2,
  cex.main = 1.4
)

dev.off()

sample_cor <- cor(
  t(datExpr),
  use = "pairwise.complete.obs"
)

sample_conn <- rowSums(
  sample_cor,
  na.rm = TRUE
)

sample_conn_z <- as.numeric(
  scale(sample_conn)
)

sample_qc <- data.frame(
  sample = rownames(datExpr),
  connectivity = sample_conn,
  connectivity_Z = sample_conn_z,
  stringsAsFactors = FALSE
)

sample_qc <- sample_qc[
  order(sample_qc$connectivity_Z),
  ,
  drop = FALSE
]

write.table(
  sample_qc,
  "TCGA_OSCC_Nature8_WGCNA_sample_QC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  datExpr,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

saveRDS(
  gene_mad,
  "TCGA_OSCC_Nature8_WGCNA_gene_MAD.rds"
)

saveRDS(
  tumor_meta,
  "TCGA_OSCC_Nature8_WGCNA_tumor_metadata.rds"
)

cat("\n========================================\n")

cat(
  "Samples:",
  nrow(datExpr),
  "\n"
)

cat(
  "Genes:",
  ncol(datExpr),
  "\n"
)

cat(
  "goodSamplesGenes allOK:",
  gsg$allOK,
  "\n"
)

cat(
  "Samples with connectivity Z < -2.5:",
  sum(
    sample_qc$connectivity_Z < -2.5
  ),
  "\n"
)

cat(
  "Samples with connectivity Z < -3:",
  sum(
    sample_qc$connectivity_Z < -3
  ),
  "\n"
)

cat("========================================\n")
