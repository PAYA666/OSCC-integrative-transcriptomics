rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

required_pkgs <- c(
  "ggplot2",
  "pheatmap"
)

pkg_check <- sapply(
  required_pkgs,
  requireNamespace,
  quietly = TRUE
)

if (!all(pkg_check)) {

  stop(
    paste0(
      "Missing package(s): ",
      paste(
        names(pkg_check)[!pkg_check],
        collapse = ", "
      )
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})

project_dir <- file.path(project_root, "Part3_Cisplatin")

step1_file <- file.path(
  project_dir,
  "Step1_FINAL_V4_OFFLINE",
  "Part3_Step1_FINAL_V4_AllResults.rds"
)

step2_file <- file.path(
  project_dir,
  "Step2_Hallmark_GSEA_V2",
  "Part3_Step2_Hallmark_V2_COMPLETE.rds"
)

if (!file.exists(step1_file)) {
  stop(
    paste0(
      "Step 1 file not found:\n",
      step1_file
    )
  )
}

if (!file.exists(step2_file)) {
  stop(
    paste0(
      "Step 2 file not found:\n",
      step2_file
    )
  )
}

out_dir <- file.path(
  project_dir,
  "Step3_Gene_Convergence"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

step1 <- readRDS(
  step1_file
)

step2 <- readRDS(
  step2_file
)

stopifnot(
  all(
    c(
      "GSE197561",
      "GSE111585"
    ) %in%
      names(step1)
  )
)

d197 <- step1$GSE197561$deg_all

d111 <- step1$GSE111585$deg_all

needed_cols <- c(
  "Gene",
  "logFC",
  "t",
  "P.Value",
  "adj.P.Val"
)

stopifnot(
  all(
    needed_cols %in%
      colnames(d197)
  )
)

stopifnot(
  all(
    needed_cols %in%
      colnames(d111)
  )
)

x197 <- d197[
  ,
  needed_cols
]

colnames(x197) <- c(
  "Gene",
  "logFC_197",
  "t_197",
  "P_197",
  "FDR_197"
)

x111 <- d111[
  ,
  needed_cols
]

colnames(x111) <- c(
  "Gene",
  "logFC_111",
  "t_111",
  "P_111",
  "FDR_111"
)

gene_meta <- merge(
  x197,
  x111,
  by = "Gene",
  all = FALSE
)

cat(
  "\nCommon genes measured in both datasets:",
  nrow(gene_meta),
  "\n"
)

gene_meta$Same_direction <-
  sign(
    gene_meta$logFC_197
  ) ==
  sign(
    gene_meta$logFC_111
  )

gene_meta$Direction <- ifelse(

  gene_meta$Same_direction &
    gene_meta$logFC_197 > 0,

  "Up_in_Resistant",

  ifelse(

    gene_meta$Same_direction &
      gene_meta$logFC_197 < 0,

    "Down_in_Resistant",

    "Opposite"
  )
)

gene_meta$Strict_197 <-
  gene_meta$FDR_197 < 0.05 &
  abs(
    gene_meta$logFC_197
  ) > 1

gene_meta$Strict_111 <-
  gene_meta$FDR_111 < 0.05 &
  abs(
    gene_meta$logFC_111
  ) > 1

gene_meta$N_Strict_DEG <-
  as.integer(
    gene_meta$Strict_197
  ) +
  as.integer(
    gene_meta$Strict_111
  )

signed_z <- function(
    p,
    effect
) {

  p2 <- pmin(
    pmax(
      p,
      1e-300
    ),
    1
  )

  z <- qnorm(
    p2 / 2,
    lower.tail = FALSE
  )

  z <- z *
    sign(effect)

  return(z)
}

gene_meta$Z_197 <- signed_z(
  gene_meta$P_197,
  gene_meta$logFC_197
)

gene_meta$Z_111 <- signed_z(
  gene_meta$P_111,
  gene_meta$logFC_111
)

gene_meta$Meta_Z <- (
  gene_meta$Z_197 +
    gene_meta$Z_111
) / sqrt(2)

gene_meta$Meta_P <- 2 * pnorm(
  abs(
    gene_meta$Meta_Z
  ),
  lower.tail = FALSE
)

gene_meta$Meta_FDR <- p.adjust(
  gene_meta$Meta_P,
  method = "BH"
)

gene_meta$Meta_Direction <- ifelse(
  gene_meta$Meta_Z > 0,
  "Up_in_Resistant",
  "Down_in_Resistant"
)

gene_meta$Mean_logFC <- rowMeans(
  gene_meta[
    ,
    c(
      "logFC_197",
      "logFC_111"
    )
  ]
)

gene_meta$Min_abs_logFC <- pmin(
  abs(
    gene_meta$logFC_197
  ),
  abs(
    gene_meta$logFC_111
  )
)

gene_meta$Evidence_Class <-
  "Opposite_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction
] <- "Same_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction &
    gene_meta$Meta_FDR < 0.05
] <- "Meta_FDR005_same_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction &
    gene_meta$N_Strict_DEG == 1
] <- "One_strict_DEG_same_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction &
    gene_meta$N_Strict_DEG == 2
] <- "Strict_DEG_both_same_direction"

gene_meta$Evidence_Class <- "Opposite_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction
] <- "Same_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction &
    gene_meta$Meta_FDR < 0.05
] <- "Meta_FDR005_same_direction"

gene_meta$Evidence_Class[
  gene_meta$Same_direction &
    gene_meta$Meta_FDR < 0.05 &
    gene_meta$N_Strict_DEG == 1
] <- "Meta_FDR005_one_strict"

gene_meta$Evidence_Class[
  gene_meta$Same_direction &
    gene_meta$N_Strict_DEG == 2
] <- "Strict_DEG_both_same_direction"

get_leading_edge <- function(
    gsea_object,
    pathway
) {

  idx <- match(
    pathway,
    gsea_object$all$pathway
  )

  if (is.na(idx)) {
    return(
      character(0)
    )
  }

  x <- gsea_object$all$leadingEdge[[idx]]

  if (is.null(x)) {
    return(
      character(0)
    )
  }

  return(
    unique(
      as.character(x)
    )
  )
}

core_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING"
)

cat("\n========== PATHWAY CHECK ==========\n")

for (p in core_pathways) {

  cat(
    p,
    "| GSE197561 =",
    p %in% step2$GSE197561$all$pathway,
    "| GSE111585 =",
    p %in% step2$GSE111585$all$pathway,
    "\n"
  )
}

leading_list <- list()

for (pathway in core_pathways) {

  le197 <- get_leading_edge(
    step2$GSE197561,
    pathway
  )

  le111 <- get_leading_edge(
    step2$GSE111585,
    pathway
  )

  all_genes <- union(
    le197,
    le111
  )

  cat(
    pathway,
    "| LE197 =",
    length(le197),
    "| LE111 =",
    length(le111),
    "| overlap =",
    length(
      intersect(
        le197,
        le111
      )
    ),
    "\n"
  )

  if (length(all_genes) == 0) {
    next
  }

  tmp <- data.frame(
    Gene = all_genes,
    Pathway = pathway,
    LeadingEdge_197 =
      all_genes %in%
      le197,
    LeadingEdge_111 =
      all_genes %in%
      le111,
    stringsAsFactors = FALSE
  )

  tmp$LeadingEdge_Both <-
    tmp$LeadingEdge_197 &
    tmp$LeadingEdge_111

  leading_list[[pathway]] <- tmp
}

if (length(leading_list) == 0) {
  stop(
    "No leading-edge genes were extracted."
  )
}

leading_long <- do.call(
  rbind,
  leading_list
)

rownames(
  leading_long
) <- NULL

write.table(
  leading_long,
  file = file.path(
    out_dir,
    "Part3_CorePathway_LeadingEdge_Long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

leading_summary <- do.call(
  rbind,
  lapply(
    core_pathways,
    function(pathway) {

      x <- leading_long[
        leading_long$Pathway ==
          pathway,
        ,
        drop = FALSE
      ]

      data.frame(
        Pathway = pathway,
        LeadingEdge_197 =
          sum(
            x$LeadingEdge_197
          ),
        LeadingEdge_111 =
          sum(
            x$LeadingEdge_111
          ),
        LeadingEdge_Union =
          nrow(x),
        LeadingEdge_Intersection =
          sum(
            x$LeadingEdge_Both
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  leading_summary,
  file = file.path(
    out_dir,
    "Part3_CorePathway_LeadingEdge_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

core_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING"
)

pathway_short_names <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFA_NFKB",
  "HALLMARK_HYPOXIA" = "Hypoxia",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_TGF_BETA_SIGNALING" = "TGFb"
)

leading_genes <- unique(leading_long$Gene)

pathway_gene_evidence <- data.frame(
  Gene = leading_genes,
  stringsAsFactors = FALSE
)

for (pathway in core_pathways) {

  short_name <- unname(pathway_short_names[pathway])

  x <- leading_long[
    leading_long$Pathway == pathway,
    ,
    drop = FALSE
  ]

  col197 <- paste0(short_name, "_LE197")
  col111 <- paste0(short_name, "_LE111")

  pathway_gene_evidence[, col197] <-
    pathway_gene_evidence$Gene %in%
    x$Gene[x$LeadingEdge_197]

  pathway_gene_evidence[, col111] <-
    pathway_gene_evidence$Gene %in%
    x$Gene[x$LeadingEdge_111]
}

le_cols <- grep(
  "_LE197$|_LE111$",
  colnames(pathway_gene_evidence),
  value = TRUE
)

pathway_gene_evidence$LeadingEdge_Evidence_Count <-
  rowSums(
    pathway_gene_evidence[
      ,
      le_cols,
      drop = FALSE
    ]
  )

pathway_gene_evidence$BothDataset_Pathway_Count <- 0L

for (short_name in c(
  "TNFA_NFKB",
  "Hypoxia",
  "EMT",
  "TGFb"
)) {

  col197 <- paste0(short_name, "_LE197")
  col111 <- paste0(short_name, "_LE111")

  both_here <-
    pathway_gene_evidence[, col197] &
    pathway_gene_evidence[, col111]

  pathway_gene_evidence$BothDataset_Pathway_Count <-
    pathway_gene_evidence$BothDataset_Pathway_Count +
    as.integer(both_here)
}

pathway_gene_evidence$Pathway_Tags <- ""

for (i in seq_len(nrow(pathway_gene_evidence))) {

  tags <- character(0)

  if (
    pathway_gene_evidence$TNFA_NFKB_LE197[i] ||
    pathway_gene_evidence$TNFA_NFKB_LE111[i]
  ) {
    tags <- c(tags, "TNFA_NFKB")
  }

  if (
    pathway_gene_evidence$Hypoxia_LE197[i] ||
    pathway_gene_evidence$Hypoxia_LE111[i]
  ) {
    tags <- c(tags, "Hypoxia")
  }

  if (
    pathway_gene_evidence$EMT_LE197[i] ||
    pathway_gene_evidence$EMT_LE111[i]
  ) {
    tags <- c(tags, "EMT")
  }

  if (
    pathway_gene_evidence$TGFb_LE197[i] ||
    pathway_gene_evidence$TGFb_LE111[i]
  ) {
    tags <- c(tags, "TGFb")
  }

  pathway_gene_evidence$Pathway_Tags[i] <-
    paste(tags, collapse = ";")
}

gene_meta$Core_Pathway_Hit <-
  gene_meta$Gene %in% pathway_gene_evidence$Gene

gene_evidence <- merge(
  gene_meta,
  pathway_gene_evidence,
  by = "Gene",
  all.x = TRUE
)

logical_cols <- grep(
  "_LE197$|_LE111$",
  colnames(gene_evidence),
  value = TRUE
)

for (nm in logical_cols) {

  gene_evidence[
    is.na(gene_evidence[, nm]),
    nm
  ] <- FALSE
}

gene_evidence$LeadingEdge_Evidence_Count[
  is.na(gene_evidence$LeadingEdge_Evidence_Count)
] <- 0

gene_evidence$BothDataset_Pathway_Count[
  is.na(gene_evidence$BothDataset_Pathway_Count)
] <- 0

gene_evidence$Pathway_Tags[
  is.na(gene_evidence$Pathway_Tags)
] <- ""

gene_evidence$Priority_Pool <-
  gene_evidence$Same_direction &
  gene_evidence$Meta_FDR < 0.05 &
  gene_evidence$N_Strict_DEG >= 1 &
  gene_evidence$Core_Pathway_Hit

gene_evidence$Ultra_Strict <-
  gene_evidence$Same_direction &
  gene_evidence$N_Strict_DEG == 2 &
  gene_evidence$Core_Pathway_Hit

gene_evidence <- gene_evidence[
  order(
    -as.integer(gene_evidence$Ultra_Strict),
    -as.integer(gene_evidence$Priority_Pool),
    -gene_evidence$BothDataset_Pathway_Count,
    -gene_evidence$LeadingEdge_Evidence_Count,
    gene_evidence$Meta_FDR,
    -abs(gene_evidence$Meta_Z)
  ),
  ,
  drop = FALSE
]

priority_pool <- gene_evidence[
  gene_evidence$Priority_Pool,
  ,
  drop = FALSE
]

ultra_strict <- gene_evidence[
  gene_evidence$Ultra_Strict,
  ,
  drop = FALSE
]

strict_both_all <- gene_evidence[
  gene_evidence$Same_direction &
  gene_evidence$N_Strict_DEG == 2,
  ,
  drop = FALSE
]

meta_concordant <- gene_evidence[
  gene_evidence$Same_direction &
  gene_evidence$Meta_FDR < 0.05,
  ,
  drop = FALSE
]

spearman_logFC <- cor(
  gene_meta$logFC_197,
  gene_meta$logFC_111,
  method = "spearman",
  use = "complete.obs"
)

spearman_t <- cor(
  gene_meta$t_197,
  gene_meta$t_111,
  method = "spearman",
  use = "complete.obs"
)

same_direction_pct <- 100 *
  mean(gene_meta$Same_direction)

concordance_summary <- data.frame(
  Common_genes = nrow(gene_meta),
  Same_direction_genes = sum(gene_meta$Same_direction),
  Same_direction_percent = same_direction_pct,
  Spearman_logFC = spearman_logFC,
  Spearman_t = spearman_t,
  Strict_DEG_both_same = nrow(strict_both_all),
  Meta_FDR005_same = nrow(meta_concordant),
  Core_leading_edge_genes = nrow(pathway_gene_evidence),
  Priority_pool = nrow(priority_pool),
  Ultra_strict = nrow(ultra_strict),
  stringsAsFactors = FALSE
)

write.table(
  concordance_summary,
  file = file.path(
    out_dir,
    "Part3_Gene_Concordance_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pathway_gene_evidence,
  file = file.path(
    out_dir,
    "Part3_CorePathway_GeneEvidence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_evidence,
  file = file.path(
    out_dir,
    "Part3_Gene_Evidence_Master.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  priority_pool,
  file = file.path(
    out_dir,
    "Part3_Priority_Gene_Pool.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ultra_strict,
  file = file.path(
    out_dir,
    "Part3_Ultra_Strict_CorePathway_Genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  strict_both_all,
  file = file.path(
    out_dir,
    "Part3_Strict_DEG_Both_SameDirection.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  meta_concordant,
  file = file.path(
    out_dir,
    "Part3_MetaFDR005_SameDirection.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

step3_results <- list(
  concordance_summary = concordance_summary,
  gene_meta = gene_meta,
  gene_evidence = gene_evidence,
  leading_edge_summary = leading_summary,
  leading_edge_long = leading_long,
  pathway_gene_evidence = pathway_gene_evidence,
  priority_pool = priority_pool,
  ultra_strict = ultra_strict,
  strict_both_all = strict_both_all,
  meta_concordant = meta_concordant
)

saveRDS(
  step3_results,
  file = file.path(
    out_dir,
    "Part3_Step3_GeneConvergence_COMPLETE.rds"
  )
)

cat("\n")
cat("####################################################\n")
cat("STEP 3 FINISHED SUCCESSFULLY\n")
cat("####################################################\n\n")

print(concordance_summary)
