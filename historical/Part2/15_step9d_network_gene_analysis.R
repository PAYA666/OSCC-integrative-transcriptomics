project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

library(WGCNA)

options(stringsAsFactors = FALSE)

expr_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

stopifnot(file.exists(expr_file))

datExpr <- readRDS(expr_file)

me_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes.rds"
)

stopifnot(file.exists(me_file))

MEs <- readRDS(me_file)

module_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv"
)

stopifnot(file.exists(module_file))

gene_module <- read.delim(
  module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

smoking_file <- file.path(
  proj_dir,
  "Part2_Smoking_sample_metadata_patient341.rds"
)

stopifnot(file.exists(smoking_file))

smoking_meta <- readRDS(smoking_file)

deg_file <- file.path(
  proj_dir,
  "Part2_Smoking_DEG_Annotated.rds"
)

stopifnot(file.exists(deg_file))

deg <- readRDS(deg_file)

wgcna_sample <- rownames(datExpr)

wgcna_patient <- substr(
  wgcna_sample,
  1,
  12
)

smoking_aligned <- smoking_meta[
  match(
    wgcna_patient,
    smoking_meta$patient_id
  ),
]

stopifnot(
  all(
    smoking_aligned$patient_id ==
      wgcna_patient
  )
)

smoking_numeric <- rep(
  NA_real_,
  nrow(smoking_aligned)
)

smoking_numeric[
  smoking_aligned$smoking_bin == "Never"
] <- 0

smoking_numeric[
  smoking_aligned$smoking_bin == "Ever"
] <- 1

complete_idx <- !is.na(
  smoking_numeric
)

datExpr_complete <- datExpr[
  complete_idx,
  ,
  drop = FALSE
]

smoking_numeric_complete <- smoking_numeric[
  complete_idx
]

if (!is.null(rownames(MEs))) {

  me_patient <- substr(
    rownames(MEs),
    1,
    12
  )

  MEs <- MEs[
    match(
      wgcna_patient,
      me_patient
    ),
    ,
    drop = FALSE
  ]
}

MEs_complete <- MEs[
  complete_idx,
  ,
  drop = FALSE
]

MM_cor <- cor(
  datExpr_complete,
  MEs_complete,
  use = "pairwise.complete.obs",
  method = "pearson"
)

MM_p <- corPvalueStudent(
  MM_cor,
  nSamples = nrow(datExpr_complete)
)

GS_cor <- cor(
  datExpr_complete,
  smoking_numeric_complete,
  use = "pairwise.complete.obs",
  method = "pearson"
)

GS_p <- corPvalueStudent(
  GS_cor,
  nSamples = nrow(datExpr_complete)
)

wgcna_gene_id <- colnames(datExpr_complete)

idx_module <- match(
  wgcna_gene_id,
  gene_module$gene_id
)

idx_deg <- match(
  wgcna_gene_id,
  deg$gene_id
)

network_gene_table <- data.frame(

  gene_id = wgcna_gene_id,

  gene_name =
    gene_module$gene_name[idx_module],

  module_numeric =
    gene_module$module_numeric[idx_module],

  module_color =
    gene_module$module_color[idx_module],

  Smoking_GS =
    as.numeric(GS_cor),

  Smoking_GS_P =
    as.numeric(GS_p),

  log2FoldChange =
    deg$log2FoldChange[idx_deg],

  DEG_padj =
    deg$padj[idx_deg],

  stringsAsFactors = FALSE
)

network_gene_table$MM <- NA_real_

network_gene_table$MM_P <- NA_real_

for (i in seq_len(nrow(network_gene_table))) {

  mod_num <- network_gene_table$module_numeric[i]

  me_name <- paste0(
    "ME",
    mod_num
  )

  if (me_name %in% colnames(MM_cor)) {

    network_gene_table$MM[i] <-
      MM_cor[
        i,
        me_name
      ]

    network_gene_table$MM_P[i] <-
      MM_p[
        i,
        me_name
      ]
  }
}

network_gene_table$Smoking_DEG <- (
  !is.na(network_gene_table$DEG_padj) &
  network_gene_table$DEG_padj < 0.05
)

network_gene_table$Strong_Smoking_DEG <- (
  network_gene_table$Smoking_DEG &
  abs(network_gene_table$log2FoldChange) >= 1
)

target_modules <- c(
  "turquoise",
  "black",
  "red"
)

target_network <- network_gene_table[
  network_gene_table$module_color %in%
    target_modules,
]

target_deg <- target_network[
  target_network$Smoking_DEG,
]

candidate_core <- target_deg[
  !is.na(target_deg$MM) &
  abs(target_deg$MM) >= 0.70,
]

candidate_core$NetworkScore <-
  abs(candidate_core$MM) *
  abs(candidate_core$Smoking_GS) *
  (-log10(
    pmax(
      candidate_core$DEG_padj,
      1e-300
    )
  ))

candidate_core <- candidate_core[
  order(
    candidate_core$NetworkScore,
    decreasing = TRUE
  ),
]

saveRDS(
  network_gene_table,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_GeneLevel_MM_GS.rds"
  )
)

write.table(
  network_gene_table,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_GeneLevel_MM_GS.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  target_deg,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_TargetModules_DEGs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  candidate_core,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_CoreCandidateGenes.rds"
  )
)

write.table(
  candidate_core,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_CoreCandidateGenes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 9D finished\n")
cat("========================================\n")

cat(
  "Target-module genes:",
  nrow(target_network),
  "\n"
)

cat(
  "Smoking-DEGs in target modules:",
  nrow(target_deg),
  "\n"
)

cat(
  "Core candidates |MM|>=0.70:",
  nrow(candidate_core),
  "\n"
)
