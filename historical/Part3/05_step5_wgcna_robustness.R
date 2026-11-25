rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

required_pkgs <- c(
  "WGCNA"
)

pkg_ok <- sapply(
  required_pkgs,
  requireNamespace,
  quietly = TRUE
)

print(pkg_ok)

if (!all(pkg_ok)) {
  stop(
    "WGCNA package is required for Step 5."
  )
}

suppressPackageStartupMessages({
  library(WGCNA)
})

root <- project_root

part3_dir <- file.path(
  root,
  "Part3_Cisplatin"
)

step3_file <- file.path(
  part3_dir,
  "Step3_Gene_Convergence",
  "Part3_Step3_GeneConvergence_COMPLETE.rds"
)

step4_file <- file.path(
  part3_dir,
  "Step4_TCGA_WGCNA_Projection",
  "Part3_Step4_TCGA_WGCNA_COMPLETE.rds"
)

datExpr_file <- file.path(
  root,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

ME_file <- file.path(
  root,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes.rds"
)

required_files <- c(
  step3_file,
  step4_file,
  datExpr_file,
  ME_file
)

if (!all(file.exists(required_files))) {

  cat("\nMissing files:\n")

  print(
    required_files[
      !file.exists(required_files)
    ]
  )

  stop(
    "Required Step 5 file missing."
  )
}

out_dir <- file.path(
  part3_dir,
  "Step5_TCGA_WGCNA_Robustness"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

step3 <- readRDS(
  step3_file
)

step4 <- readRDS(
  step4_file
)

datExpr <- readRDS(
  datExpr_file
)

MEs <- readRDS(
  ME_file
)

expr_symbol <- as.matrix(
  step4$expression_symbol_341
)

module_symbol <- as.data.frame(
  step4$module_symbol_table
)

priority_table <- as.data.frame(
  step3$priority_pool
)

strict_table <- as.data.frame(
  step3$strict_both_all
)

ultra_table <- as.data.frame(
  step3$ultra_strict
)

sample_ids <- colnames(
  expr_symbol
)

stopifnot(
  length(sample_ids) == 341
)

stopifnot(
  all(
    sample_ids %in%
      rownames(MEs)
  )
)

stopifnot(
  all(
    sample_ids %in%
      rownames(datExpr)
  )
)

MEs341 <- as.data.frame(
  MEs[
    sample_ids,
    ,
    drop = FALSE
  ]
)

datExpr341 <- as.data.frame(
  datExpr[
    sample_ids,
    ,
    drop = FALSE
  ]
)

module_color_map <- unique(
  module_symbol[
    ,
    c(
      "module_numeric",
      "module_color"
    ),
    drop = FALSE
  ]
)

module_color_map$module_numeric <-
  as.integer(
    module_color_map$module_numeric
  )

module_color_map <- module_color_map[
  order(
    module_color_map$module_numeric
  ),
  ,
  drop = FALSE
]

me_map <- data.frame(
  ME = colnames(MEs341),
  module_numeric = as.integer(
    sub(
      "^ME",
      "",
      colnames(MEs341)
    )
  ),
  stringsAsFactors = FALSE
)

me_map$module_color <-
  module_color_map$module_color[
    match(
      me_map$module_numeric,
      module_color_map$module_numeric
    )
  ]

cat("\n========== ME MAP ==========\n")
print(me_map)

make_score <- function(
    gene_table,
    genes_keep,
    expression_matrix
) {

  genes_keep <- unique(
    as.character(
      genes_keep
    )
  )

  genes_keep <- genes_keep[
    genes_keep %in%
      rownames(expression_matrix)
  ]

  if (length(genes_keep) < 3) {

    return(
      rep(
        NA_real_,
        ncol(expression_matrix)
      )
    )
  }

  gene_info <- gene_table[
    match(
      genes_keep,
      gene_table$Gene
    ),
    ,
    drop = FALSE
  ]

  x <- expression_matrix[
    genes_keep,
    ,
    drop = FALSE
  ]

  xz <- t(
    scale(
      t(x)
    )
  )

  good <- apply(
    xz,
    1,
    function(z) {
      all(
        is.finite(z)
      )
    }
  )

  xz <- xz[
    good,
    ,
    drop = FALSE
  ]

  if (nrow(xz) < 3) {

    return(
      rep(
        NA_real_,
        ncol(expression_matrix)
      )
    )
  }

  gene_info <- gene_info[
    match(
      rownames(xz),
      gene_info$Gene
    ),
    ,
    drop = FALSE
  ]

  direction_sign <- ifelse(
    gene_info$Meta_Z >= 0,
    1,
    -1
  )

  signed_x <- sweep(
    xz,
    1,
    direction_sign,
    FUN = "*"
  )

  score <- colMeans(
    signed_x,
    na.rm = TRUE
  )

  names(score) <- colnames(
    expression_matrix
  )

  return(score)
}

priority_genes <- unique(
  as.character(
    priority_table$Gene
  )
)

priority_genes <- priority_genes[
  priority_genes %in%
    rownames(expr_symbol)
]

strict_genes <- unique(
  as.character(
    strict_table$Gene
  )
)

strict_genes <- strict_genes[
  strict_genes %in%
    rownames(expr_symbol)
]

cat(
  "\nPriority genes available in TCGA:",
  length(priority_genes),
  "\n"
)

cat(
  "StrictBoth genes available in TCGA:",
  length(strict_genes),
  "\n"
)

priority_original <- make_score(
  gene_table = priority_table,
  genes_keep = priority_genes,
  expression_matrix = expr_symbol
)

strict_original <- make_score(
  gene_table = strict_table,
  genes_keep = strict_genes,
  expression_matrix = expr_symbol
)

run_module_depletion <- function(
    gene_table,
    original_genes,
    score_name
) {

  output <- data.frame()

  for (i in seq_len(nrow(me_map))) {

    me_name <- me_map$ME[i]

    module_color <- me_map$module_color[i]

    module_genes <- module_symbol$gene_name[
      module_symbol$module_color ==
        module_color
    ]

    removed_genes <- intersect(
      original_genes,
      module_genes
    )

    remaining_genes <- setdiff(
      original_genes,
      removed_genes
    )

    depleted_score <- make_score(
      gene_table = gene_table,
      genes_keep = remaining_genes,
      expression_matrix = expr_symbol
    )

    me_vector <- MEs341[
      sample_ids,
      me_name
    ]

    original_score <- make_score(
      gene_table = gene_table,
      genes_keep = original_genes,
      expression_matrix = expr_symbol
    )

    ct_original <- suppressWarnings(
      cor.test(
        original_score[sample_ids],
        me_vector,
        method = "spearman",
        exact = FALSE
      )
    )

    ct_depleted <- suppressWarnings(
      cor.test(
        depleted_score[sample_ids],
        me_vector,
        method = "spearman",
        exact = FALSE
      )
    )

    tmp <- data.frame(
      Score = score_name,
      ME = me_name,
      Module = module_color,
      Original_N = length(original_genes),
      Removed_N = length(removed_genes),
      Remaining_N = length(remaining_genes),
      Original_Rho = unname(
        ct_original$estimate
      ),
      Original_P = ct_original$p.value,
      Depleted_Rho = unname(
        ct_depleted$estimate
      ),
      Depleted_P = ct_depleted$p.value,
      stringsAsFactors = FALSE
    )

    tmp$Delta_Abs_Rho <-
      abs(tmp$Depleted_Rho) -
      abs(tmp$Original_Rho)

    tmp$Retained_Fraction <-
      ifelse(
        abs(tmp$Original_Rho) > 0,
        abs(tmp$Depleted_Rho) /
          abs(tmp$Original_Rho),
        NA_real_
      )

    tmp$Removed_Genes <-
      paste(
        removed_genes,
        collapse = ";"
      )

    output <- rbind(
      output,
      tmp
    )
  }

  output$Depleted_FDR <- p.adjust(
    output$Depleted_P,
    method = "BH"
  )

  output <- output[
    order(
      output$Depleted_FDR,
      -abs(
        output$Depleted_Rho
      )
    ),
    ,
    drop = FALSE
  ]

  return(output)
}

priority_depletion <- run_module_depletion(
  gene_table = priority_table,
  original_genes = priority_genes,
  score_name = "Priority52"
)

strict_depletion <- run_module_depletion(
  gene_table = strict_table,
  original_genes = strict_genes,
  score_name = "StrictBoth97"
)

module_depletion_all <- rbind(
  priority_depletion,
  strict_depletion
)

write.table(
  module_depletion_all,
  file = file.path(
    out_dir,
    "Part3_ModuleDepleted_Score_ME_Association.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

priority_map <- step4$candidate_module_mapping

priority_map <- priority_map[
  priority_map$Gene_Set ==
    "Priority52",
  ,
  drop = FALSE
]

priority_map <- priority_map[
  !duplicated(
    priority_map$Gene
  ),
  ,
  drop = FALSE
]

candidate_network <- merge(
  priority_table,
  priority_map[
    ,
    c(
      "Gene",
      "In_TCGA_VST",
      "In_WGCNA_top5000",
      "WGCNA_gene_id",
      "module_numeric",
      "module_color"
    ),
    drop = FALSE
  ],
  by = "Gene",
  all.x = TRUE
)

candidate_network$Assigned_ME <-
  ifelse(
    !is.na(
      candidate_network$module_numeric
    ),
    paste0(
      "ME",
      candidate_network$module_numeric
    ),
    NA_character_
  )

candidate_network$Assigned_kME <- NA_real_

for (i in seq_len(nrow(candidate_network))) {

  gene <- candidate_network$Gene[i]

  me_name <- candidate_network$Assigned_ME[i]

  ensembl_id <- candidate_network$WGCNA_gene_id[i]

  if (
    is.na(me_name) ||
    is.na(ensembl_id) ||
    !(me_name %in% colnames(MEs341)) ||
    !(ensembl_id %in% colnames(datExpr341))
  ) {
    next
  }

  gene_vector <- datExpr341[
    sample_ids,
    ensembl_id
  ]

  me_vector <- MEs341[
    sample_ids,
    me_name
  ]

  candidate_network$Assigned_kME[i] <-
    as.numeric(
      WGCNA::bicor(
        gene_vector,
        me_vector,
        use = "pairwise.complete.obs"
      )
    )
}

candidate_network$LOO_Rho <- NA_real_
candidate_network$LOO_P <- NA_real_

for (i in seq_len(nrow(candidate_network))) {

  gene <- candidate_network$Gene[i]

  if (
    !(gene %in% rownames(expr_symbol))
  ) {
    next
  }

  genes_without_gene <- setdiff(
    priority_genes,
    gene
  )

  loo_score <- make_score(
    gene_table = priority_table,
    genes_keep = genes_without_gene,
    expression_matrix = expr_symbol
  )

  gene_vector <- expr_symbol[
    gene,
    sample_ids
  ]

  ct <- suppressWarnings(
    cor.test(
      as.numeric(gene_vector),
      loo_score[sample_ids],
      method = "spearman",
      exact = FALSE
    )
  )

  candidate_network$LOO_Rho[i] <-
    unname(
      ct$estimate
    )

  candidate_network$LOO_P[i] <-
    ct$p.value
}

candidate_network$LOO_FDR <- p.adjust(
  candidate_network$LOO_P,
  method = "BH"
)

candidate_network$Ultra_Strict <-
  candidate_network$Gene %in%
  ultra_table$Gene

candidate_network$Network_Interpretable <-
  !is.na(
    candidate_network$module_color
  ) &
  candidate_network$module_color !=
    "grey"

candidate_network <- candidate_network[
  order(
    -as.integer(
      candidate_network$Ultra_Strict
    ),
    -as.integer(
      candidate_network$Network_Interpretable
    ),
    -abs(
      candidate_network$Assigned_kME
    ),
    candidate_network$Meta_FDR,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

write.table(
  candidate_network,
  file = file.path(
    out_dir,
    "Part3_Priority52_TCGA_WGCNA_CandidateEvidence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ultra_network <- candidate_network[
  candidate_network$Ultra_Strict,
  ,
  drop = FALSE
]

write.table(
  ultra_network,
  file = file.path(
    out_dir,
    "Part3_Ultra7_TCGA_WGCNA_Robustness.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

priority_non_grey <- priority_depletion[
  priority_depletion$Module != "grey",
  ,
  drop = FALSE
]

strongest_priority_module <-
  priority_non_grey$Module[
    which.max(
      abs(
        priority_non_grey$Depleted_Rho
      )
    )
  ]

strongest_priority_rho <-
  priority_non_grey$Depleted_Rho[
    which.max(
      abs(
        priority_non_grey$Depleted_Rho
      )
    )
  ]

summary_step5 <- data.frame(
  TCGA_samples = length(sample_ids),
  Priority_genes_used = length(priority_genes),
  Strict_genes_used = length(strict_genes),
  Priority_modules_depleted_tested =
    nrow(priority_depletion),
  Strongest_nonGrey_module_after_depletion =
    strongest_priority_module,
  Strongest_nonGrey_depleted_rho =
    strongest_priority_rho,
  Priority_candidates_with_WGCNA =
    sum(
      candidate_network$In_WGCNA_top5000,
      na.rm = TRUE
    ),
  Ultra_candidates =
    nrow(ultra_network),
  Ultra_nonGrey =
    sum(
      ultra_network$Network_Interpretable,
      na.rm = TRUE
    ),
  stringsAsFactors = FALSE
)

write.table(
  summary_step5,
  file = file.path(
    out_dir,
    "Part3_Step5_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

step5_results <- list(
  summary = summary_step5,
  me_map = me_map,
  priority_depletion = priority_depletion,
  strict_depletion = strict_depletion,
  candidate_network = candidate_network,
  ultra_network = ultra_network
)

saveRDS(
  step5_results,
  file = file.path(
    out_dir,
    "Part3_Step5_TCGA_WGCNA_Robustness_COMPLETE.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step5_sessionInfo.txt"
  )
)

cat("\n")
cat("####################################################\n")
cat("STEP 5 FINISHED SUCCESSFULLY\n")
cat("####################################################\n\n")

print(summary_step5)
