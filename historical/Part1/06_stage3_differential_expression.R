rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(DESeq2)
library(dplyr)
library(ggplot2)

count_df <- readRDS(
  "TCGA_OSCC_Nature8_raw_count_matrix.rds"
)

meta <- readRDS(
  "TCGA_OSCC_Nature8_metadata.rds"
)

cat("Raw matrix:", dim(count_df), "\n")
cat("Metadata:", dim(meta), "\n")

print(table(meta$sample_type))

count_pc <- count_df[
  count_df$gene_type == "protein_coding",
  ,
  drop = FALSE
]

cat("\nProtein-coding genes:", nrow(count_pc), "\n")

sample_ids <- meta$sample_barcode

stopifnot(
  all(sample_ids %in% colnames(count_pc))
)

count_mat <- as.matrix(
  count_pc[, sample_ids, drop = FALSE]
)

storage.mode(count_mat) <- "integer"

rownames(count_mat) <- count_pc$gene_id

gene_anno <- count_pc[, c(
  "gene_id",
  "gene_name",
  "gene_type"
)]

rownames(gene_anno) <- gene_anno$gene_id

keep <- rowSums(count_mat >= 10) >= 10

count_filtered <- count_mat[
  keep,
  ,
  drop = FALSE
]

gene_anno_filtered <- gene_anno[
  rownames(count_filtered),
  ,
  drop = FALSE
]

cat("\nGenes before filtering:",
    nrow(count_mat), "\n")

cat("Genes after filtering:",
    nrow(count_filtered), "\n")

coldata <- data.frame(
  sample_barcode = meta$sample_barcode,
  sample_type = meta$sample_type,
  row.names = meta$sample_barcode,
  stringsAsFactors = FALSE
)

coldata$sample_type <- factor(
  coldata$sample_type,
  levels = c(
    "Solid Tissue Normal",
    "Primary Tumor"
  )
)

stopifnot(
  identical(
    colnames(count_filtered),
    rownames(coldata)
  )
)

dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData = coldata,
  design = ~ sample_type
)

dds <- DESeq(
  dds
)

res <- results(
  dds,
  contrast = c(
    "sample_type",
    "Primary Tumor",
    "Solid Tissue Normal"
  ),
  alpha = 0.05
)

res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

res_df <- res_df %>%
  left_join(
    gene_anno_filtered,
    by = "gene_id"
  )

res_df <- res_df %>%
  select(
    gene_id,
    gene_name,
    gene_type,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  )

res_df$DEG_status <- "NS"

res_df$DEG_status[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange > 1
] <- "Up"

res_df$DEG_status[
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange < -1
] <- "Down"

deg_df <- res_df %>%
  filter(
    DEG_status %in% c(
      "Up",
      "Down"
    )
  )

deg_up <- deg_df %>%
  filter(DEG_status == "Up")

deg_down <- deg_df %>%
  filter(DEG_status == "Down")

res_df <- res_df %>%
  arrange(padj)

deg_df <- deg_df %>%
  arrange(padj)

vsd <- vst(
  dds,
  blind = FALSE
)

vst_mat <- assay(vsd)

volcano_df <- res_df

volcano_df$minus_log10_padj <-
  -log10(volcano_df$padj)

finite_max <- max(
  volcano_df$minus_log10_padj[
    is.finite(volcano_df$minus_log10_padj)
  ],
  na.rm = TRUE
)

volcano_df$minus_log10_padj[
  is.infinite(volcano_df$minus_log10_padj)
] <- finite_max + 1

p_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = DEG_status
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.2
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Down" = "#377EB8",
      "NS" = "grey70",
      "Up" = "#E41A1C"
    )
  ) +
  labs(
    title = "TCGA Nature8 OSCC",
    subtitle = "Primary Tumor vs Solid Tissue Normal",
    x = "log2 Fold Change",
    y = "-log10 adjusted P value",
    color = NULL
  ) +
  theme_classic(base_size = 13)

ggsave(
  "Figure_Nature8_OSCC_DEG_Volcano.png",
  p_volcano,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  "Figure_Nature8_OSCC_DEG_Volcano.pdf",
  p_volcano,
  width = 7,
  height = 6
)

write.table(
  res_df,
  "TCGA_OSCC_Nature8_DESeq2_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  deg_df,
  "TCGA_OSCC_Nature8_DEG.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  deg_up,
  "TCGA_OSCC_Nature8_DEG_UP.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  deg_down,
  "TCGA_OSCC_Nature8_DEG_DOWN.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  dds,
  "TCGA_OSCC_Nature8_DESeq2_dds.rds"
)

saveRDS(
  res_df,
  "TCGA_OSCC_Nature8_DESeq2_results.rds"
)

saveRDS(
  deg_df,
  "TCGA_OSCC_Nature8_DEG_results.rds"
)

saveRDS(
  vsd,
  "TCGA_OSCC_Nature8_VST.rds"
)

saveRDS(
  vst_mat,
  "TCGA_OSCC_Nature8_VST_matrix.rds"
)

cat("\n========================================\n")

cat(
  "Tumor samples:",
  sum(coldata$sample_type == "Primary Tumor"),
  "\n"
)

cat(
  "Normal samples:",
  sum(coldata$sample_type == "Solid Tissue Normal"),
  "\n"
)

cat(
  "Protein-coding genes before filtering:",
  nrow(count_mat),
  "\n"
)

cat(
  "Genes entering DESeq2:",
  nrow(count_filtered),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(deg_df),
  "\n"
)

cat(
  "Up:",
  nrow(deg_up),
  "\n"
)

cat(
  "Down:",
  nrow(deg_down),
  "\n"
)

cat("========================================\n")
cat("Stage 3 completed.\n")
