rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

required_pkgs <- c(
  "AnnotationDbi",
  "org.Hs.eg.db",
  "limma"
)

pkg_ok <- sapply(
  required_pkgs,
  requireNamespace,
  quietly = TRUE
)

print(pkg_ok)

if (!all(pkg_ok)) {

  stop(
    paste0(
      "Missing package(s): ",
      paste(
        names(pkg_ok)[!pkg_ok],
        collapse = ", "
      )
    )
  )
}

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(limma)
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

datExpr_file <- file.path(
  root,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

ME_file <- file.path(
  root,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes.rds"
)

module_file <- file.path(
  root,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv"
)

vst_file <- file.path(
  root,
  "TCGA_OSCC_Nature8_VST_matrix.rds"
)

all_files <- c(
  step3_file,
  datExpr_file,
  ME_file,
  module_file,
  vst_file
)

if (!all(file.exists(all_files))) {

  cat("\nMissing files:\n")

  print(
    all_files[
      !file.exists(all_files)
    ]
  )

  stop(
    "Required file missing."
  )
}

out_dir <- file.path(
  part3_dir,
  "Step4_TCGA_WGCNA_Projection"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

step3 <- readRDS(
  step3_file
)

datExpr <- readRDS(
  datExpr_file
)

MEs <- readRDS(
  ME_file
)

module_table <- read.delim(
  module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

vst <- readRDS(
  vst_file
)

stopifnot(
  nrow(datExpr) == 341
)

stopifnot(
  nrow(MEs) == 341
)

stopifnot(
  ncol(datExpr) == 5000
)

stopifnot(
  all(
    c(
      "gene_id",
      "module_numeric",
      "module_color",
      "gene_name"
    ) %in%
      colnames(module_table)
  )
)

sample_ids <- rownames(
  datExpr
)

if (
  !all(
    sample_ids %in%
      colnames(vst)
  )
) {

  stop(
    "Not all WGCNA samples are present in VST matrix."
  )
}

if (
  !all(
    sample_ids %in%
      rownames(MEs)
  )
) {

  stop(
    "Not all WGCNA samples are present in module eigengenes."
  )
}

vst341 <- as.matrix(
  vst[
    ,
    sample_ids,
    drop = FALSE
  ]
)

MEs341 <- as.data.frame(
  MEs[
    sample_ids,
    ,
    drop = FALSE
  ]
)

stopifnot(
  identical(
    colnames(vst341),
    sample_ids
  )
)

stopifnot(
  identical(
    rownames(MEs341),
    sample_ids
  )
)

cat(
  "\nTCGA tumor samples aligned:",
  length(sample_ids),
  "\n"
)

vst_ids <- rownames(
  vst341
)

vst_ensembl <- sub(
  "\\..*$",
  "",
  vst_ids
)

ensembl_fraction <- mean(
  grepl(
    "^ENSG",
    vst_ensembl
  )
)

cat(
  "Fraction of VST rows looking like Ensembl IDs:",
  round(
    ensembl_fraction * 100,
    2
  ),
  "%\n"
)

if (
  ensembl_fraction < 0.50
) {

  stop(
    "VST row names do not appear to be Ensembl IDs."
  )
}

unique_ensembl <- unique(
  vst_ensembl
)

symbol_map <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique_ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
)

vst_symbols <- unname(
  symbol_map[
    vst_ensembl
  ]
)

mapped_vst <- !is.na(
  vst_symbols
) &
  vst_symbols != ""

cat(
  "VST Ensembl -> Symbol mapped:",
  sum(mapped_vst),
  "/",
  length(mapped_vst),
  "\n"
)

expr_symbol <- limma::avereps(
  vst341[
    mapped_vst,
    ,
    drop = FALSE
  ],
  ID = vst_symbols[
    mapped_vst
  ]
)

expr_symbol <- as.matrix(
  expr_symbol
)

gene_sd <- apply(
  expr_symbol,
  1,
  sd,
  na.rm = TRUE
)

expr_symbol <- expr_symbol[
  is.finite(gene_sd) &
    gene_sd > 0,
  ,
  drop = FALSE
]

cat(
  "TCGA symbol-level expression matrix:",
  nrow(expr_symbol),
  "genes x",
  ncol(expr_symbol),
  "tumors\n"
)

priority_table <- step3$priority_pool

ultra_table <- step3$ultra_strict

strict_table <- step3$strict_both_all

cat(
  "\nExternal gene sets:\n"
)

cat(
  "Priority =",
  nrow(priority_table),
  "\n"
)

cat(
  "Ultra-strict =",
  nrow(ultra_table),
  "\n"
)

cat(
  "Strict-both =",
  nrow(strict_table),
  "\n"
)

make_signed_score <- function(
    gene_table,
    score_name,
    expression_matrix
) {

  genes <- unique(
    as.character(
      gene_table$Gene
    )
  )

  available <- genes[
    genes %in%
      rownames(
        expression_matrix
      )
  ]

  missing <- setdiff(
    genes,
    available
  )

  cat(
    "\n",
    score_name,
    ": mapped ",
    length(available),
    "/",
    length(genes),
    " genes\n",
    sep = ""
  )

  if (
    length(missing) > 0
  ) {

    cat(
      score_name,
      " missing genes: ",
      paste(
        missing,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }

  if (
    length(available) < 3
  ) {

    stop(
      paste0(
        score_name,
        " has fewer than 3 genes available in TCGA."
      )
    )
  }

  x <- expression_matrix[
    available,
    ,
    drop = FALSE
  ]

  x_z <- t(
    scale(
      t(x)
    )
  )

  good <- apply(
    x_z,
    1,
    function(z) {
      all(
        is.finite(z)
      )
    }
  )

  x_z <- x_z[
    good,
    ,
    drop = FALSE
  ]

  available <- rownames(
    x_z
  )

  meta_z <- gene_table$Meta_Z[
    match(
      available,
      gene_table$Gene
    )
  ]

  direction_sign <- ifelse(
    meta_z >= 0,
    1,
    -1
  )

  signed_z <- sweep(
    x_z,
    1,
    direction_sign,
    FUN = "*"
  )

  score <- colMeans(
    signed_z,
    na.rm = TRUE
  )

  names(score) <- colnames(
    signed_z
  )

  return(
    list(
      score = score,
      genes_used = available,
      genes_missing = missing,
      n_used = length(available)
    )
  )
}

priority_score_result <- make_signed_score(
  gene_table = priority_table,
  score_name = "Priority52",
  expression_matrix = expr_symbol
)

strict_score_result <- make_signed_score(
  gene_table = strict_table,
  score_name = "StrictBoth97",
  expression_matrix = expr_symbol
)

ultra_score_result <- make_signed_score(
  gene_table = ultra_table,
  score_name = "Ultra7",
  expression_matrix = expr_symbol
)

score_table <- data.frame(
  Sample = sample_ids,
  Priority52_Score =
    priority_score_result$score[
      sample_ids
    ],
  StrictBoth97_Score =
    strict_score_result$score[
      sample_ids
    ],
  Ultra7_Score =
    ultra_score_result$score[
      sample_ids
    ],
  stringsAsFactors = FALSE
)

rownames(
  score_table
) <- score_table$Sample

write.table(
  score_table,
  file = file.path(
    out_dir,
    "Part3_TCGA_ResistanceAssociated_Scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

score_cor <- cor(
  score_table[
    ,
    c(
      "Priority52_Score",
      "StrictBoth97_Score",
      "Ultra7_Score"
    )
  ],
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.table(
  score_cor,
  file = file.path(
    out_dir,
    "Part3_TCGA_Score_Spearman_Correlation.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

module_color_map <- unique(
  module_table[
    ,
    c(
      "module_numeric",
      "module_color"
    )
  ]
)

module_color_map <- module_color_map[
  order(
    module_color_map$module_numeric
  ),
  ,
  drop = FALSE
]

score_me_cor <- data.frame()

score_names <- c(
  "Priority52_Score",
  "StrictBoth97_Score",
  "Ultra7_Score"
)

for (score_name in score_names) {

  score_vector <- score_table[
    sample_ids,
    score_name
  ]

  for (me_name in colnames(MEs341)) {

    me_vector <- MEs341[
      ,
      me_name
    ]

    ct <- suppressWarnings(
      cor.test(
        score_vector,
        me_vector,
        method = "spearman",
        exact = FALSE
      )
    )

    tmp <- data.frame(
      Score = score_name,
      ME = me_name,
      Rho = unname(
        ct$estimate
      ),
      P_value = ct$p.value,
      stringsAsFactors = FALSE
    )

    score_me_cor <- rbind(
      score_me_cor,
      tmp
    )
  }
}

score_me_cor$FDR <- NA_real_

for (
  score_name in unique(
    score_me_cor$Score
  )
) {

  ii <- score_me_cor$Score ==
    score_name

  score_me_cor$FDR[ii] <- p.adjust(
    score_me_cor$P_value[ii],
    method = "BH"
  )
}

score_me_cor$module_numeric <- as.integer(
  sub(
    "^ME",
    "",
    score_me_cor$ME
  )
)

score_me_cor$module_color <-
  module_color_map$module_color[
    match(
      score_me_cor$module_numeric,
      module_color_map$module_numeric
    )
  ]

score_me_cor <- score_me_cor[
  order(
    score_me_cor$Score,
    score_me_cor$FDR,
    -abs(
      score_me_cor$Rho
    )
  ),
  ,
  drop = FALSE
]

write.table(
  score_me_cor,
  file = file.path(
    out_dir,
    "Part3_TCGA_Score_WGCNA_ME_Association.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_mad <- apply(
  datExpr,
  2,
  mad,
  na.rm = TRUE
)

module_symbol <- module_table

module_symbol$gene_name <- trimws(
  as.character(
    module_symbol$gene_name
  )
)

module_symbol <- module_symbol[
  !is.na(
    module_symbol$gene_name
  ) &
    module_symbol$gene_name != "",
  ,
  drop = FALSE
]

module_symbol$WGCNA_MAD <-
  gene_mad[
    match(
      module_symbol$gene_id,
      names(gene_mad)
    )
  ]

mad_for_sort <- module_symbol$WGCNA_MAD

mad_for_sort[
  is.na(
    mad_for_sort
  )
] <- -Inf

module_symbol <- module_symbol[
  order(
    module_symbol$gene_name,
    -mad_for_sort
  ),
  ,
  drop = FALSE
]

duplicate_symbol_count <- sum(
  duplicated(
    module_symbol$gene_name
  )
)

cat(
  "\nDuplicated gene symbols in WGCNA table:",
  duplicate_symbol_count,
  "\n"
)

module_symbol_unique <- module_symbol[
  !duplicated(
    module_symbol$gene_name
  ),
  ,
  drop = FALSE
]

cat(
  "Unique gene symbols in WGCNA universe:",
  nrow(module_symbol_unique),
  "\n"
)

map_candidate_set <- function(
    gene_table,
    set_name
) {

  genes <- unique(
    as.character(
      gene_table$Gene
    )
  )

  idx <- match(
    genes,
    module_symbol_unique$gene_name
  )

  out <- data.frame(
    Gene = genes,
    Gene_Set = set_name,
    In_TCGA_VST =
      genes %in%
      rownames(expr_symbol),
    In_WGCNA_top5000 =
      !is.na(idx),
    WGCNA_gene_id =
      module_symbol_unique$gene_id[idx],
    module_numeric =
      module_symbol_unique$module_numeric[idx],
    module_color =
      module_symbol_unique$module_color[idx],
    stringsAsFactors = FALSE
  )

  return(out)
}

priority_module_map <- map_candidate_set(
  priority_table,
  "Priority52"
)

ultra_module_map <- map_candidate_set(
  ultra_table,
  "Ultra7"
)

strict_module_map <- map_candidate_set(
  strict_table,
  "StrictBoth97"
)

candidate_module_map <- rbind(
  priority_module_map,
  ultra_module_map,
  strict_module_map
)

write.table(
  candidate_module_map,
  file = file.path(
    out_dir,
    "Part3_Candidate_to_WGCNA_Module_Mapping.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

module_enrichment <- function(
    genes,
    gene_set_name
) {

  genes <- unique(
    as.character(
      genes
    )
  )

  universe <- module_symbol_unique

  universe$Hit <-
    universe$gene_name %in%
    genes

  colors <- sort(
    unique(
      universe$module_color
    )
  )

  result <- data.frame()

  for (color in colors) {

    in_module <-
      universe$module_color ==
      color

    a <- sum(
      in_module &
        universe$Hit
    )

    b <- sum(
      in_module &
        !universe$Hit
    )

    c <- sum(
      !in_module &
        universe$Hit
    )

    d <- sum(
      !in_module &
        !universe$Hit
    )

    tab <- matrix(
      c(
        a,
        b,
        c,
        d
      ),
      nrow = 2,
      byrow = TRUE
    )

    ft <- fisher.test(
      tab,
      alternative = "greater"
    )

    tmp <- data.frame(
      Gene_Set = gene_set_name,
      Module = color,
      Module_Size = sum(in_module),
      GeneSet_in_WGCNA = sum(
        universe$Hit
      ),
      Hits_in_Module = a,
      Odds_Ratio = unname(
        ft$estimate
      ),
      P_value = ft$p.value,
      stringsAsFactors = FALSE
    )

    result <- rbind(
      result,
      tmp
    )
  }

  result$FDR <- p.adjust(
    result$P_value,
    method = "BH"
  )

  result <- result[
    order(
      result$FDR,
      -result$Odds_Ratio
    ),
    ,
    drop = FALSE
  ]

  return(result)
}

enrich_strict <- module_enrichment(
  genes = strict_table$Gene,
  gene_set_name = "StrictBoth97"
)

enrich_priority <- module_enrichment(
  genes = priority_table$Gene,
  gene_set_name = "Priority52"
)

module_enrichment_all <- rbind(
  enrich_strict,
  enrich_priority
)

write.table(
  module_enrichment_all,
  file = file.path(
    out_dir,
    "Part3_ResistanceGenes_WGCNA_Module_Enrichment.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_step4 <- data.frame(
  TCGA_tumors = length(sample_ids),

  VST_symbol_genes =
    nrow(expr_symbol),

  Priority_total =
    nrow(priority_table),

  Priority_in_TCGA =
    priority_score_result$n_used,

  StrictBoth_total =
    nrow(strict_table),

  StrictBoth_in_TCGA =
    strict_score_result$n_used,

  Ultra_total =
    nrow(ultra_table),

  Ultra_in_TCGA =
    ultra_score_result$n_used,

  Priority_in_WGCNA =
    sum(
      priority_module_map$In_WGCNA_top5000
    ),

  StrictBoth_in_WGCNA =
    sum(
      strict_module_map$In_WGCNA_top5000
    ),

  Ultra_in_WGCNA =
    sum(
      ultra_module_map$In_WGCNA_top5000
    ),

  stringsAsFactors = FALSE
)

write.table(
  summary_step4,
  file = file.path(
    out_dir,
    "Part3_Step4_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

step4_results <- list(
  summary = summary_step4,
  expression_symbol_341 = expr_symbol,
  scores = score_table,
  score_correlations = score_cor,
  score_ME_association = score_me_cor,
  module_symbol_table = module_symbol_unique,
  candidate_module_mapping = candidate_module_map,
  module_enrichment = module_enrichment_all,
  priority_score_genes =
    priority_score_result$genes_used,
  strict_score_genes =
    strict_score_result$genes_used,
  ultra_score_genes =
    ultra_score_result$genes_used
)

saveRDS(
  step4_results,
  file = file.path(
    out_dir,
    "Part3_Step4_TCGA_WGCNA_COMPLETE.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step4_sessionInfo.txt"
  )
)

cat("\n")
cat("####################################################\n")
cat("STEP 4 FINISHED SUCCESSFULLY\n")
cat("####################################################\n\n")

print(
  summary_step4
)
