rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(WGCNA)

options(stringsAsFactors = FALSE)

allowWGCNAThreads()

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

cat("Input dimensions:\n")
print(dim(datExpr))

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

stopifnot(gsg$allOK)

softPower <- 8

set.seed(12345)

net <- blockwiseModules(
  datExpr,

  power = softPower,

  corType = "bicor",

  networkType = "signed",
  TOMType = "signed",

  minModuleSize = 30,

  deepSplit = 2,

  mergeCutHeight = 0.25,

  pamRespectsDendro = FALSE,

  numericLabels = TRUE,

  maxBlockSize = 6000,

  saveTOMs = TRUE,

  saveTOMFileBase =
    "TCGA_OSCC_Nature8_WGCNA_TOM",

  verbose = 3
)

moduleColors <- labels2colors(
  net$colors
)

cat("\nModule number including grey:\n")
print(length(unique(moduleColors)))

cat("\nModule sizes:\n")
print(
  sort(
    table(moduleColors),
    decreasing = TRUE
  )
)

MEs <- orderMEs(
  net$MEs
)

cat("\nModule eigengene matrix:\n")
print(dim(MEs))

pdf(
  "Figure_Nature8_OSCC_WGCNA_ModuleDendrogram.pdf",
  width = 14,
  height = 8
)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[
    net$blockGenes[[1]]
  ],
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main =
    "WGCNA gene dendrogram and module assignment"
)

dev.off()

png(
  "Figure_Nature8_OSCC_WGCNA_ModuleDendrogram.png",
  width = 2800,
  height = 1600,
  res = 200
)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[
    net$blockGenes[[1]]
  ],
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main =
    "WGCNA gene dendrogram and module assignment"
)

dev.off()

gene_module <- data.frame(

  gene_id = colnames(datExpr),

  module_numeric = net$colors,

  module_color = moduleColors,

  stringsAsFactors = FALSE
)

res_df <- readRDS(
  "TCGA_OSCC_Nature8_DESeq2_results.rds"
)

gene_annotation <- unique(
  res_df[
    ,
    c(
      "gene_id",
      "gene_name"
    )
  ]
)

gene_module <- merge(
  gene_module,
  gene_annotation,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

gene_module <- gene_module[
  match(
    colnames(datExpr),
    gene_module$gene_id
  ),
  ,
  drop = FALSE
]

module_size <- as.data.frame(
  table(moduleColors)
)

colnames(module_size) <- c(
  "module",
  "gene_count"
)

module_size <- module_size[
  order(
    module_size$gene_count,
    decreasing = TRUE
  ),
]

saveRDS(
  net,
  "TCGA_OSCC_Nature8_WGCNA_network_power8.rds"
)

saveRDS(
  moduleColors,
  "TCGA_OSCC_Nature8_WGCNA_moduleColors.rds"
)

saveRDS(
  MEs,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes.rds"
)

saveRDS(
  gene_module,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.rds"
)

write.table(
  gene_module,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  module_size,
  "TCGA_OSCC_Nature8_WGCNA_module_size.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
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
  "Power:",
  softPower,
  "\n"
)

cat(
  "Modules including grey:",
  length(unique(moduleColors)),
  "\n"
)

cat(
  "Genes assigned to grey:",
  sum(moduleColors == "grey"),
  "\n"
)

cat("========================================\n")
