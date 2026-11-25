project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

key_genes <- c(

  "CYP1A1",
  "CES1",

  "SFRP4",
  "SFRP2",
  "PDGFRA",
  "OGN",
  "FNDC1",
  "TK1",
  "MYBL2"
)

sample_info <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_sample_metadata_patient341.rds"
  )
)

stopifnot(
  nrow(sample_info) == 341,
  length(unique(sample_info$patient_id)) == 341
)

sample_info333 <- sample_info[
  !is.na(sample_info$smoking_bin),
]

stopifnot(
  nrow(sample_info333) == 333
)

cat("Smoking-complete unique patients:\n")
print(table(sample_info333$smoking_bin))

expr_vst <- readRDS(
  file.path(
    proj_dir,
    "TCGA_OSCC_Nature8_VST_matrix.rds"
  )
)

missing_samples <- setdiff(
  sample_info333$sample_barcode,
  colnames(expr_vst)
)

if (length(missing_samples) > 0) {
  stop(
    paste(
      "Missing VST samples:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

expr333 <- expr_vst[
  ,
  sample_info333$sample_barcode,
  drop = FALSE
]

stopifnot(
  ncol(expr333) == 333
)

stopifnot(
  identical(
    colnames(expr333),
    sample_info333$sample_barcode
  )
)

counts_raw <- readRDS(
  file.path(
    proj_dir,
    "TCGA_OSCC_Nature8_raw_count_matrix.rds"
  )
)

gene_annotation <- unique(
  counts_raw[
    ,
    c(
      "gene_id",
      "gene_name",
      "gene_type"
    )
  ]
)

idx_anno <- match(
  rownames(expr333),
  gene_annotation$gene_id
)

vst_symbols <- gene_annotation$gene_name[
  idx_anno
]

names(vst_symbols) <- rownames(expr333)

key_gene_ids <- sapply(
  key_genes,
  function(g) {

    hits <- names(vst_symbols)[
      !is.na(vst_symbols) &
      vst_symbols == g
    ]

    if (length(hits) == 0) {
      return(NA_character_)
    }

    hits[1]
  }
)

if (any(is.na(key_gene_ids))) {

  cat(
    "Genes not found:\n"
  )

  print(
    names(key_gene_ids)[
      is.na(key_gene_ids)
    ]
  )
}

key_gene_ids <- key_gene_ids[
  !is.na(key_gene_ids)
]

key_expr <- expr333[
  key_gene_ids,
  ,
  drop = FALSE
]

rownames(key_expr) <- names(key_gene_ids)

Never_idx <- (
  sample_info333$smoking_bin == "Never"
)

Ever_idx <- (
  sample_info333$smoking_bin == "Ever"
)

expression_summary <- data.frame(

  gene_name = rownames(key_expr),

  Never_mean = rowMeans(
    key_expr[, Never_idx, drop = FALSE]
  ),

  Ever_mean = rowMeans(
    key_expr[, Ever_idx, drop = FALSE]
  ),

  Never_median = apply(
    key_expr[, Never_idx, drop = FALSE],
    1,
    median
  ),

  Ever_median = apply(
    key_expr[, Ever_idx, drop = FALSE],
    1,
    median
  ),

  stringsAsFactors = FALSE
)

deg <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_Annotated.rds"
  )
)

formal_stats <- deg[
  deg$gene_name %in% key_genes,
  c(
    "gene_name",
    "gene_type",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "DEG_status"
  )
]

formal_stats$gene_name <- factor(
  formal_stats$gene_name,
  levels = key_genes
)

formal_stats <- formal_stats[
  order(formal_stats$gene_name),
]

formal_stats$gene_name <- as.character(
  formal_stats$gene_name
)

key_gene_summary <- merge(
  formal_stats,
  expression_summary,
  by = "gene_name",
  all.x = TRUE,
  sort = FALSE
)

key_gene_summary$gene_name <- factor(
  key_gene_summary$gene_name,
  levels = key_genes
)

key_gene_summary <- key_gene_summary[
  order(key_gene_summary$gene_name),
]

key_gene_summary$gene_name <- as.character(
  key_gene_summary$gene_name
)

key_gene_summary$Role <- ifelse(
  key_gene_summary$gene_name %in%
    c("CYP1A1", "CES1"),
  "Exposure_response",
  "Network_candidate"
)

expr_long <- do.call(
  rbind,
  lapply(
    rownames(key_expr),
    function(g) {

      data.frame(
        gene_name = g,
        sample_barcode = colnames(key_expr),
        patient_id = sample_info333$patient_id,
        smoking_bin = sample_info333$smoking_bin,
        VST_expression = as.numeric(
          key_expr[g,]
        ),
        stringsAsFactors = FALSE
      )

    }
  )
)

write.table(
  key_gene_summary,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_DESeq2_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  key_gene_summary,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_DESeq2_FINAL.rds"
  )
)

write.table(
  expr_long,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_VST_LongFormat.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  expr_long,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_VST_LongFormat.rds"
  )
)

cat("\n========================================\n")
cat("Step 12A FINAL finished\n")
cat("========================================\n")

cat(
  "Patients:",
  nrow(sample_info333),
  "\n"
)

cat(
  "Key genes:",
  nrow(key_gene_summary),
  "\n"
)

print(
  key_gene_summary[
    ,
    c(
      "gene_name",
      "Role",
      "log2FoldChange",
      "padj",
      "Never_mean",
      "Ever_mean"
    )
  ]
)
