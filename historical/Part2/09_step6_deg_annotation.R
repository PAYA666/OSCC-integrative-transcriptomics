project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

library(DESeq2)

count_matrix <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_count_matrix_patient341.rds"
  )
)

sample_info <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_sample_metadata_patient341.rds"
  )
)

keep <- !is.na(
  sample_info$smoking_bin
)

count_matrix_deg <-
  count_matrix[,keep]

sample_info_deg <-
  sample_info[keep,]

rownames(sample_info_deg) <-
  sample_info_deg$sample_barcode

count_matrix_deg <- round(
  count_matrix_deg
)

dds <- DESeqDataSetFromMatrix(

  countData = count_matrix_deg,

  colData = sample_info_deg,

  design = ~ smoking_bin

)

keep_gene <-
  rowSums(
    counts(dds) >= 10
  ) >=
  floor(
    ncol(dds)*0.1
  )

dds <- dds[keep_gene,]

cat(
"Genes after filtering:",
nrow(dds),
"\n"
)

dds <- DESeq(dds)

res <- results(
  dds,
  contrast=c(
    "smoking_bin",
    "Ever",
    "Never"
  )
)

res <- res[
  order(res$padj),
]

deg <- as.data.frame(res)

saveRDS(
  deg,
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_DESeq2_Ever_vs_Never.rds"
  )
)

write.table(
  deg,
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_DESeq2_Ever_vs_Never.tsv"
  ),
  sep="\t",
  quote=FALSE
)

cat("DESeq2 finished\n")
