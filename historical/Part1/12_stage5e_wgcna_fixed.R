rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(WGCNA)
library(survival)
library(dplyr)

options(stringsAsFactors = FALSE)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

moduleColors <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_moduleColors.rds"
)

MEs_color <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes_COLOR.rds"
)

gene_module <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.rds"
)

clinical_wgcna <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_clinical_aligned.rds"
)

deg_df <- readRDS(
  "TCGA_OSCC_Nature8_DEG_results.rds"
)

candidate_modules <- c(
  "blue",
  "black"
)

candidate_gene_ids <- gene_module$gene_id[
  gene_module$module_color %in% candidate_modules
]

cat("\nCandidate module gene numbers:\n")
print(
  table(
    gene_module$module_color[
      gene_module$module_color %in% candidate_modules
    ]
  )
)

MM_table <- data.frame(
  gene_id = colnames(datExpr),
  module_color = moduleColors,
  stringsAsFactors = FALSE
)

for (mod in candidate_modules) {

  me_name <- paste0("ME", mod)

  mm <- bicor(
    datExpr,
    MEs_color[, me_name],
    use = "pairwise.complete.obs"
  )

  MM_table[[paste0("MM_", mod)]] <- as.numeric(mm)
}

anno <- unique(
  gene_module[
    ,
    c(
      "gene_id",
      "gene_name"
    )
  ]
)

MM_table <- dplyr::left_join(
  MM_table,
  anno,
  by = "gene_id"
)

surv_keep <- (
  !is.na(clinical_wgcna$OS_time) &
  clinical_wgcna$OS_time > 0 &
  !is.na(clinical_wgcna$OS_event)
)

surv_dat <- clinical_wgcna[
  surv_keep,
  ,
  drop = FALSE
]

expr_surv <- datExpr[
  surv_keep,
  ,
  drop = FALSE
]

cat("\nPatients entering gene-level Cox:\n")
print(nrow(expr_surv))

cox_gene_list <- vector(
  "list",
  length(candidate_gene_ids)
)

names(cox_gene_list) <- candidate_gene_ids

for (i in seq_along(candidate_gene_ids)) {

  gene <- candidate_gene_ids[i]

  tmp <- data.frame(
    OS_time = surv_dat$OS_time,
    OS_event = surv_dat$OS_event,
    expr = expr_surv[, gene]
  )

  fit <- try(
    survival::coxph(
      survival::Surv(
        OS_time,
        OS_event
      ) ~ expr,
      data = tmp
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {

    cox_gene_list[[i]] <- data.frame(
      gene_id = gene,
      beta = NA_real_,
      HR = NA_real_,
      lower95 = NA_real_,
      upper95 = NA_real_,
      P = NA_real_
    )

  } else {

    s <- summary(fit)

    cox_gene_list[[i]] <- data.frame(
      gene_id = gene,
      beta = s$coefficients[1, "coef"],
      HR = s$coefficients[1, "exp(coef)"],
      lower95 = s$conf.int[1, "lower .95"],
      upper95 = s$conf.int[1, "upper .95"],
      P = s$coefficients[1, "Pr(>|z|)"]
    )
  }
}

gene_cox <- dplyr::bind_rows(
  cox_gene_list
)

gene_cox$FDR <- p.adjust(
  gene_cox$P,
  method = "BH"
)

gene_cox <- dplyr::left_join(
  gene_cox,
  MM_table,
  by = "gene_id"
)

deg_info <- dplyr::select(
  deg_df,
  gene_id,
  log2FoldChange,
  padj,
  DEG_status
)

gene_cox <- dplyr::left_join(
  gene_cox,
  deg_info,
  by = "gene_id"
)

gene_cox$MM <- ifelse(
  gene_cox$module_color == "blue",
  gene_cox$MM_blue,
  ifelse(
    gene_cox$module_color == "black",
    gene_cox$MM_black,
    NA_real_
  )
)

candidate_final <- gene_cox %>%
  dplyr::filter(
    abs(MM) >= 0.70,
    P < 0.05,
    DEG_status %in% c("Up", "Down")
  ) %>%
  dplyr::arrange(P)

candidate_fdr <- gene_cox %>%
  dplyr::filter(
    abs(MM) >= 0.70,
    FDR < 0.05,
    DEG_status %in% c("Up", "Down")
  ) %>%
  dplyr::arrange(FDR)

write.table(
  gene_cox,
  "TCGA_OSCC_Nature8_WGCNA_candidateModules_geneCox_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  candidate_final,
  "TCGA_OSCC_Nature8_WGCNA_DEG_Cox_MM_candidates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  candidate_fdr,
  "TCGA_OSCC_Nature8_WGCNA_DEG_Cox_MM_FDR_candidates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  gene_cox,
  "TCGA_OSCC_Nature8_WGCNA_candidateModules_geneCox_ALL.rds"
)

saveRDS(
  candidate_final,
  "TCGA_OSCC_Nature8_WGCNA_DEG_Cox_MM_candidates.rds"
)

saveRDS(
  candidate_fdr,
  "TCGA_OSCC_Nature8_WGCNA_DEG_Cox_MM_FDR_candidates.rds"
)

cat("\n========================================\n")

cat(
  "Blue genes:",
  sum(gene_cox$module_color == "blue"),
  "\n"
)

cat(
  "Black genes:",
  sum(gene_cox$module_color == "black"),
  "\n"
)

cat(
  "Gene Cox P < 0.05:",
  sum(gene_cox$P < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Gene Cox FDR < 0.05:",
  sum(gene_cox$FDR < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Final exploratory candidates:",
  nrow(candidate_final),
  "\n"
)

cat(
  "Strict FDR candidates:",
  nrow(candidate_fdr),
  "\n"
)

cat("========================================\n")
