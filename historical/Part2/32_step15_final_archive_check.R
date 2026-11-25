project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root
fig_dir  <- file.path(
  proj_dir,
  "Figures_Final"
)

exposure <- readRDS(
  file.path(
    proj_dir,
    "Part2_TCGA_OSCC_Nature8_Exposure_341_CLEAN.rds"
  )
)

deg <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_Annotated.rds"
  )
)

gsea <- as.data.frame(
  readRDS(
    file.path(
      proj_dir,
      "Part2_Smoking_Hallmark_GSEA_MAIN.rds"
    )
  )
)

keygenes <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_DESeq2_FINAL.rds"
  )
)

status <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_StatusPattern_FINAL.rds"
  )
)

fib <- read.delim(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_FibroblastCorrelation.tsv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

N_total <- nrow(exposure)

N_never <- sum(
  exposure$smoking_bin == "Never",
  na.rm = TRUE
)

N_ever <- sum(
  exposure$smoking_bin == "Ever",
  na.rm = TRUE
)

N_missing <- sum(
  is.na(exposure$smoking_bin)
)

N_DEG <- sum(
  deg$padj < 0.05,
  na.rm = TRUE
)

N_strong_up <- sum(
  deg$padj < 0.05 &
  deg$log2FoldChange >= 1,
  na.rm = TRUE
)

N_strong_down <- sum(
  deg$padj < 0.05 &
  deg$log2FoldChange <= -1,
  na.rm = TRUE
)

N_GSEA_sig <- sum(
  gsea$p.adjust < 0.05,
  na.rm = TRUE
)

N_GSEA_positive <- sum(
  gsea$p.adjust < 0.05 &
  gsea$NES > 0,
  na.rm = TRUE
)

N_GSEA_negative <- sum(
  gsea$p.adjust < 0.05 &
  gsea$NES < 0,
  na.rm = TRUE
)

N_key_DEG <- sum(
  keygenes$padj < 0.05,
  na.rm = TRUE
)

N_current_sig <- sum(
  status$Current_vs_Never_FDR < 0.05,
  na.rm = TRUE
)

N_former_sig <- sum(
  status$Former_vs_Never_FDR < 0.05,
  na.rm = TRUE
)

N_packyear_sig <- sum(
  status$PackYears_FDR < 0.05,
  na.rm = TRUE
)

N_fibro <- sum(
  fib$Spearman_Rho > 0.50 &
  fib$FDR < 0.05,
  na.rm = TRUE
)

multi_genes <- c(
  "SFRP4",
  "PDGFRA",
  "SFRP2",
  "TK1",
  "OGN",
  "MYBL2",
  "FNDC1"
)

N_multi <- length(
  multi_genes
)

critical_genes <- c(
  "CYP1A1",
  "CES1",
  "SFRP2",
  "PDGFRA",
  "SFRP4",
  "FNDC1",
  "OGN",
  "TK1",
  "MYBL2"
)

missing_critical <- setdiff(
  critical_genes,
  keygenes$gene_name
)

figure_files <- c(

  "Figure1C_Volcano_PUBLICATION.pdf",

  "Figure2A_SmokingStatus_KeyGenes_PUBLICATION.pdf",

  "Figure2B_FibroblastCorrelation_PUBLICATION.pdf",

  "Figure3A_Smoking_Hallmark_GSEA_MAIN.pdf",

  "Figure3B_Representative_GSEA_Curves_PUBLICATION_FIXED.pdf",

  "Figure4A_WGCNA_SmokingDEG_Enrichment_PUBLICATION.pdf",

  "Figure4B_MultiEvidence_Candidates_PUBLICATION.pdf",

  "Figure5A_EvidenceConvergenceMatrix_PUBLICATION.pdf",

  "Figure6_WorkingModel_PUBLICATION_FINAL.pdf"
)

figure_exists <- file.exists(
  file.path(
    fig_dir,
    figure_files
  )
)

figure_qc <- data.frame(

  Figure = figure_files,

  Exists = figure_exists,

  stringsAsFactors = FALSE
)

qc_summary <- data.frame(

  Metric = c(

    "Total_OSCC_patients",

    "Smoking_complete",

    "Never_smokers",

    "Ever_smokers",

    "Smoking_missing",

    "Smoking_DEGs_FDR_lt_0.05",

    "Strong_up_DEGs",

    "Strong_down_DEGs",

    "Hallmark_GSEA_FDR_lt_0.05",

    "Hallmark_Ever_enriched",

    "Hallmark_Never_enriched",

    "Key_genes_DEG_FDR_lt_0.05",

    "Current_vs_Never_keygenes",

    "Former_vs_Never_keygenes",

    "PackYears_keygenes_FDR_lt_0.05",

    "Fibroblast_associated_keygenes",

    "MultiEvidence_candidates"

  ),

  Value = c(

    N_total,

    N_never + N_ever,

    N_never,

    N_ever,

    N_missing,

    N_DEG,

    N_strong_up,

    N_strong_down,

    N_GSEA_sig,

    N_GSEA_positive,

    N_GSEA_negative,

    N_key_DEG,

    N_current_sig,

    N_former_sig,

    N_packyear_sig,

    N_fibro,

    N_multi

  ),

  stringsAsFactors = FALSE
)

write.table(
  qc_summary,
  file.path(
    proj_dir,
    "FINAL_PROJECT_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  figure_qc,
  file.path(
    proj_dir,
    "FINAL_FIGURE_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

final_snapshot <- list(

  cohort = list(
    total = N_total,
    smoking_complete = N_never + N_ever,
    never = N_never,
    ever = N_ever
  ),

  DEG = list(
    significant = N_DEG,
    strong_up = N_strong_up,
    strong_down = N_strong_down
  ),

  GSEA = list(
    significant = N_GSEA_sig,
    positive = N_GSEA_positive,
    negative = N_GSEA_negative
  ),

  keygenes = keygenes,

  smoking_status = status,

  fibroblast = fib,

  multi_evidence_genes = multi_genes,

  figure_qc = figure_qc
)

saveRDS(
  final_snapshot,
  file.path(
    proj_dir,
    "FINAL_Smoking_OSCC_ResultSnapshot.rds"
  )
)

cat("\n========================================\n")
cat("FINAL PROJECT QC\n")
cat("========================================\n")

print(
  qc_summary
)

cat("\nFigure QC:\n")

print(
  figure_qc
)

cat("\nMissing critical genes:\n")

print(
  missing_critical
)

cat("\n========================================\n")
cat("FINAL RESULTS FROZEN\n")
cat("========================================\n")
