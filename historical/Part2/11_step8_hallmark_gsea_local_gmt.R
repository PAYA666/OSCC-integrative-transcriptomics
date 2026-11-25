project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

library(clusterProfiler)

gsea_rank <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_GSEA_rank_WaldStat.rds"
  )
)

gsea_rank <- sort(
  gsea_rank,
  decreasing = TRUE
)

cat(
"Genes in ranked list:",
length(gsea_rank),
"\n"
)

gmt_file <- file.path(
  proj_dir,
  "h.all.v2026.1.Hs.symbols.gmt"
)

stopifnot(
  file.exists(gmt_file)
)

hallmark <- read.gmt(
  gmt_file
)

cat(
"Hallmark geneset rows:",
nrow(hallmark),
"\n"
)

cat(
"Hallmark pathways:",
length(unique(hallmark$term)),
"\n"
)

set.seed(12345)

gsea_obj <- GSEA(

  geneList = gsea_rank,

  TERM2GENE = hallmark,

  minGSSize = 10,

  maxGSSize = 500,

  pvalueCutoff = 1,

  pAdjustMethod = "BH",

  verbose = FALSE

)

gsea_df <- as.data.frame(gsea_obj)

gsea_df$Direction <- ifelse(
  gsea_df$NES > 0,
  "Ever_enriched",
  "Never_enriched"
)

write.table(
  gsea_df,
  file.path(
    proj_dir,
    "Part2_Smoking_Hallmark_GSEA_MAIN.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

saveRDS(
  gsea_obj,
  file.path(
    proj_dir,
    "Part2_Smoking_Hallmark_GSEA_MAIN.rds"
  )
)

cat(
"Finished\n"
)
