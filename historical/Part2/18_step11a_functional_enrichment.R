project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    "clusterProfiler",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {

  BiocManager::install(
    "org.Hs.eg.db",
    ask = FALSE,
    update = FALSE
  )
}

library(clusterProfiler)
library(org.Hs.eg.db)

multi_file <- file.path(
  proj_dir,
  "Part2_Smoking_PPI_WGCNA_MultiEvidence.tsv"
)

stopifnot(file.exists(multi_file))

multi <- read.delim(
  multi_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_genes <- unique(
  multi$gene_name
)

cat(
  "Candidate genes:",
  length(candidate_genes),
  "\n"
)

print(candidate_genes)

gene_entrez <- bitr(
  candidate_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat(
  "Mapped genes:",
  nrow(gene_entrez),
  "\n"
)

ego_candidate <- enrichGO(
  gene = gene_entrez$ENTREZID,

  OrgDb = org.Hs.eg.db,

  keyType = "ENTREZID",

  ont = "BP",

  pAdjustMethod = "BH",

  pvalueCutoff = 0.05,

  qvalueCutoff = 0.2,

  readable = TRUE
)

ekk_candidate <- enrichKEGG(
  gene = gene_entrez$ENTREZID,

  organism = "hsa",

  pAdjustMethod = "BH",

  pvalueCutoff = 0.05
)

saveRDS(
  ego_candidate,
  file.path(
    proj_dir,
    "Part2_Smoking_Candidate_GO_BP.rds"
  )
)

saveRDS(
  ekk_candidate,
  file.path(
    proj_dir,
    "Part2_Smoking_Candidate_KEGG.rds"
  )
)

write.table(
  as.data.frame(ego_candidate),
  file.path(
    proj_dir,
    "Part2_Smoking_Candidate_GO_BP.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  as.data.frame(ekk_candidate),
  file.path(
    proj_dir,
    "Part2_Smoking_Candidate_KEGG.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

target_file <- file.path(
  proj_dir,
  "Part2_Smoking_WGCNA_TargetModules_DEGs.tsv"
)

target_deg <- read.delim(
  target_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

target_symbols <- unique(
  target_deg$gene_name
)

target_entrez <- bitr(
  target_symbols,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)

ego_target <- enrichGO(
  gene = target_entrez$ENTREZID,

  OrgDb = org.Hs.eg.db,

  keyType="ENTREZID",

  ont="BP",

  pAdjustMethod="BH",

  pvalueCutoff=0.05,

  qvalueCutoff=0.2,

  readable=TRUE
)

saveRDS(
  ego_target,
  file.path(
    proj_dir,
    "Part2_Smoking_TargetModule_GO_BP.rds"
  )
)

write.table(
  as.data.frame(ego_target),
  file.path(
    proj_dir,
    "Part2_Smoking_TargetModule_GO_BP.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("\n=============================\n")
cat("Step 11A finished\n")
cat("=============================\n")

cat(
  "Candidate GO terms:",
  nrow(as.data.frame(ego_candidate)),
  "\n"
)

cat(
  "Target module GO terms:",
  nrow(as.data.frame(ego_target)),
  "\n"
)
