project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

deg_file <- file.path(
  proj_dir,
  "Part2_Smoking_DEG_DESeq2_Ever_vs_Never.rds"
)

stopifnot(file.exists(deg_file))

deg <- readRDS(deg_file)

count_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_raw_count_matrix.rds"
)

stopifnot(file.exists(count_file))

counts_raw <- readRDS(count_file)

gene_anno <- counts_raw[, c(
  "gene_id",
  "gene_name",
  "gene_type"
)]

deg$gene_id <- rownames(deg)

idx <- match(
  deg$gene_id,
  gene_anno$gene_id
)

deg$gene_name <- gene_anno$gene_name[idx]
deg$gene_type <- gene_anno$gene_type[idx]

deg_annotated <- deg[, c(
  "gene_id",
  "gene_name",
  "gene_type",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

deg_annotated$DEG_status <- "NS"

deg_annotated$DEG_status[
  !is.na(deg_annotated$padj) &
  deg_annotated$padj < 0.05 &
  deg_annotated$log2FoldChange > 0
] <- "Up"

deg_annotated$DEG_status[
  !is.na(deg_annotated$padj) &
  deg_annotated$padj < 0.05 &
  deg_annotated$log2FoldChange < 0
] <- "Down"

deg_annotated$Strong_DEG <- FALSE

deg_annotated$Strong_DEG[
  !is.na(deg_annotated$padj) &
  deg_annotated$padj < 0.05 &
  abs(deg_annotated$log2FoldChange) >= 1
] <- TRUE

smoking_deg <- deg_annotated[
  !is.na(deg_annotated$padj) &
  deg_annotated$padj < 0.05,
]

smoking_strong_deg <- deg_annotated[
  deg_annotated$Strong_DEG,
]

gsea_df <- deg_annotated[
  !is.na(deg_annotated$gene_name) &
  deg_annotated$gene_name != "" &
  !is.na(deg_annotated$stat),
  c(
    "gene_name",
    "stat",
    "log2FoldChange",
    "padj"
  )
]

gsea_df$abs_stat <- abs(gsea_df$stat)

gsea_df <- gsea_df[
  order(
    gsea_df$gene_name,
    -gsea_df$abs_stat
  ),
]

gsea_df <- gsea_df[
  !duplicated(gsea_df$gene_name),
]

gsea_df <- gsea_df[
  order(
    gsea_df$stat,
    decreasing = TRUE
  ),
]

gsea_df$abs_stat <- NULL

gsea_rank <- gsea_df$stat

names(gsea_rank) <- gsea_df$gene_name

gsea_rank <- sort(
  gsea_rank,
  decreasing = TRUE
)

saveRDS(
  deg_annotated,
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_Annotated.rds"
  )
)

write.table(
  deg_annotated,
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_Annotated.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  smoking_deg,
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_FDR005.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  smoking_strong_deg,
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_FDR005_absLFC1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  gsea_rank,
  file.path(
    proj_dir,
    "Part2_Smoking_GSEA_rank_WaldStat.rds"
  )
)

write.table(
  gsea_df,
  file.path(
    proj_dir,
    "Part2_Smoking_GSEA_rank_WaldStat.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 7 finished successfully\n")
cat("========================================\n")

cat(
  "Annotated genes:",
  nrow(deg_annotated),
  "\n"
)

cat(
  "FDR < 0.05 DEGs:",
  nrow(smoking_deg),
  "\n"
)

cat(
  "Strong DEGs:",
  nrow(smoking_strong_deg),
  "\n"
)

cat(
  "Genes available for GSEA:",
  length(gsea_rank),
  "\n"
)
