rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)

res_df <- readRDS(
  "TCGA_OSCC_Nature8_DESeq2_results.rds"
)

deg_df <- readRDS(
  "TCGA_OSCC_Nature8_DEG_results.rds"
)

cat("All tested genes:", nrow(res_df), "\n")
cat("DEGs:", nrow(deg_df), "\n")
print(table(deg_df$DEG_status))

all_symbols <- unique(
  res_df$gene_name[
    !is.na(res_df$gene_name) &
      res_df$gene_name != ""
  ]
)

up_symbols <- unique(
  deg_df$gene_name[
    deg_df$DEG_status == "Up" &
      !is.na(deg_df$gene_name)
  ]
)

down_symbols <- unique(
  deg_df$gene_name[
    deg_df$DEG_status == "Down" &
      !is.na(deg_df$gene_name)
  ]
)

run_GO <- function(genes, ontology) {

  enrichGO(
    gene          = genes,
    universe      = all_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
}

go_up_bp   <- run_GO(up_symbols,   "BP")
go_up_cc   <- run_GO(up_symbols,   "CC")
go_up_mf   <- run_GO(up_symbols,   "MF")

go_down_bp <- run_GO(down_symbols, "BP")
go_down_cc <- run_GO(down_symbols, "CC")
go_down_mf <- run_GO(down_symbols, "MF")

add_info <- function(x, direction, ontology) {

  df <- as.data.frame(x)

  if (nrow(df) == 0) return(df)

  df$Direction <- direction
  df$Ontology  <- ontology

  df
}

go_list <- list(
  add_info(go_up_bp,   "Up",   "BP"),
  add_info(go_up_cc,   "Up",   "CC"),
  add_info(go_up_mf,   "Up",   "MF"),
  add_info(go_down_bp, "Down", "BP"),
  add_info(go_down_cc, "Down", "CC"),
  add_info(go_down_mf, "Down", "MF")
)

go_all <- bind_rows(go_list)

write.table(
  go_all,
  "TCGA_OSCC_Nature8_GO_ORA_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(go_all) > 0) {

  go_top <- go_all %>%
    filter(!is.na(p.adjust)) %>%
    group_by(Direction, Ontology) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 5) %>%
    ungroup()

  go_top$Description <- factor(
    go_top$Description,
    levels = rev(unique(go_top$Description))
  )

  p_go <- ggplot(
    go_top,
    aes(
      x = -log10(p.adjust),
      y = Description,
      size = Count,
      color = Direction
    )
  ) +
    geom_point(alpha = 0.85) +
    facet_grid(
      Ontology ~ Direction,
      scales = "free_y",
      space = "free_y"
    ) +
    labs(
      title = "GO enrichment of Nature8 OSCC DEGs",
      x = "-log10 adjusted P value",
      y = NULL,
      size = "Gene count",
      color = "Direction"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    "Figure_Nature8_OSCC_GO_Combined.png",
    p_go,
    width = 12,
    height = 10,
    dpi = 600
  )

  ggsave(
    "Figure_Nature8_OSCC_GO_Combined.pdf",
    p_go,
    width = 12,
    height = 10
  )
}

all_map <- bitr(
  all_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_map <- bitr(
  up_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_map <- bitr(
  down_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

all_entrez  <- unique(all_map$ENTREZID)
up_entrez   <- unique(up_map$ENTREZID)
down_entrez <- unique(down_map$ENTREZID)

kegg_up <- enrichKEGG(
  gene          = up_entrez,
  universe      = all_entrez,
  organism      = "hsa",
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20
)

kegg_down <- enrichKEGG(
  gene          = down_entrez,
  universe      = all_entrez,
  organism      = "hsa",
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.20
)

kegg_up <- setReadable(
  kegg_up,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_down <- setReadable(
  kegg_down,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_up_df <- as.data.frame(kegg_up)
kegg_down_df <- as.data.frame(kegg_down)

if (nrow(kegg_up_df) > 0)
  kegg_up_df$Direction <- "Up"

if (nrow(kegg_down_df) > 0)
  kegg_down_df$Direction <- "Down"

kegg_all <- bind_rows(
  kegg_up_df,
  kegg_down_df
)

write.table(
  kegg_all,
  "TCGA_OSCC_Nature8_KEGG_ORA_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(kegg_all) > 0) {

  kegg_top <- kegg_all %>%
    filter(!is.na(p.adjust)) %>%
    group_by(Direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()

  kegg_top$Description <- factor(
    kegg_top$Description,
    levels = rev(unique(kegg_top$Description))
  )

  p_kegg <- ggplot(
    kegg_top,
    aes(
      x = -log10(p.adjust),
      y = Description,
      size = Count,
      color = Direction
    )
  ) +
    geom_point(alpha = 0.85) +
    facet_wrap(
      ~ Direction,
      scales = "free_y"
    ) +
    labs(
      title = "KEGG enrichment of Nature8 OSCC DEGs",
      x = "-log10 adjusted P value",
      y = NULL,
      size = "Gene count",
      color = "Direction"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    "Figure_Nature8_OSCC_KEGG_Combined.png",
    p_kegg,
    width = 12,
    height = 7,
    dpi = 600
  )

  ggsave(
    "Figure_Nature8_OSCC_KEGG_Combined.pdf",
    p_kegg,
    width = 12,
    height = 7
  )
}

gsea_df <- res_df %>%
  filter(
    !is.na(stat),
    !is.na(gene_name),
    gene_name != ""
  )

gsea_df <- gsea_df %>%
  mutate(abs_stat = abs(stat)) %>%
  arrange(desc(abs_stat)) %>%
  distinct(
    gene_name,
    .keep_all = TRUE
  )

gene_list <- gsea_df$stat
names(gene_list) <- gsea_df$gene_name

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

hallmark_file <-
  file.path(project_root, "h.all.v2026.1.Hs.symbols.gmt")

stopifnot(
  file.exists(hallmark_file)
)

hallmark_gmt <- read.gmt(
  hallmark_file
)

set.seed(12345)

gsea_hallmark <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = hallmark_gmt,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = FALSE,
  seed          = TRUE
)

gsea_df_all <- as.data.frame(
  gsea_hallmark
)

if (nrow(gsea_df_all) > 0) {

  gsea_df_all$Direction <- ifelse(
    gsea_df_all$NES > 0,
    "Tumor",
    "Normal"
  )

}

gsea_sig <- gsea_df_all %>%
  filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  )

write.table(
  gsea_df_all,
  "TCGA_OSCC_Nature8_GSEA_Hallmark_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gsea_sig,
  "TCGA_OSCC_Nature8_GSEA_Hallmark_SIGNIFICANT.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  gsea_hallmark,
  "TCGA_OSCC_Nature8_GSEA_Hallmark_results.rds"
)

if (nrow(gsea_sig) > 0) {

  tumor_top <- gsea_sig %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    slice_head(n = 10)

  normal_top <- gsea_sig %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    slice_head(n = 10)

  hallmark_top <- bind_rows(
    tumor_top,
    normal_top
  )

  hallmark_top$Description <- gsub(
    "^HALLMARK_",
    "",
    hallmark_top$Description
  )

  hallmark_top$Description <- gsub(
    "_",
    " ",
    hallmark_top$Description
  )

  hallmark_top$Description <- factor(
    hallmark_top$Description,
    levels = hallmark_top$Description[
      order(hallmark_top$NES)
    ]
  )

  p_gsea <- ggplot(
    hallmark_top,
    aes(
      x = NES,
      y = Description,
      color = Direction,
      size = -log10(p.adjust)
    )
  ) +
    geom_point(alpha = 0.9) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = "Hallmark GSEA in Nature8 OSCC",
      subtitle = "Positive NES: Tumor enriched; Negative NES: Normal enriched",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      color = NULL,
      size = "-log10 FDR"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    "Figure_Nature8_OSCC_GSEA_Hallmark.png",
    p_gsea,
    width = 10,
    height = 7,
    dpi = 600
  )

  ggsave(
    "Figure_Nature8_OSCC_GSEA_Hallmark.pdf",
    p_gsea,
    width = 10,
    height = 7
  )
}

rank_table <- data.frame(
  gene_name = names(gene_list),
  Wald_stat = as.numeric(gene_list)
)

write.table(
  rank_table,
  "TCGA_OSCC_Nature8_GSEA_ranked_gene_list.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  go_all,
  "TCGA_OSCC_Nature8_GO_ORA_results.rds"
)

saveRDS(
  kegg_all,
  "TCGA_OSCC_Nature8_KEGG_ORA_results.rds"
)

cat("\n========================================\n")
cat("Stage 4 finished.\n")

cat(
  "GO significant terms:",
  nrow(go_all),
  "\n"
)

cat(
  "KEGG significant pathways:",
  nrow(kegg_all),
  "\n"
)

cat(
  "Genes in Hallmark GSEA ranked list:",
  length(gene_list),
  "\n"
)

cat(
  "Significant Hallmark pathways (FDR < 0.05):",
  nrow(gsea_sig),
  "\n"
)

cat("========================================\n")
