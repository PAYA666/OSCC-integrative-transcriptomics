rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("hgu133plus2.db", quietly = TRUE))
  BiocManager::install("hgu133plus2.db")

library(hgu133plus2.db)
library(AnnotationDbi)
library(dplyr)

expr_geo <- readRDS(
  "GSE41613_expression_local.rds"
)

pheno_geo <- readRDS(
  "GSE41613_phenotype_local.rds"
)

vital_raw <- as.character(
  pheno_geo$characteristics_ch1.4
)

vital <- sub(
  "^vital:\\s*",
  "",
  vital_raw,
  ignore.case = TRUE
)

OS_event <- ifelse(
  grepl("^Alive$", vital, ignore.case = TRUE),
  0,
  ifelse(
    grepl("^Dead", vital, ignore.case = TRUE),
    1,
    NA
  )
)

fu_raw <- as.character(
  pheno_geo$characteristics_ch1.5
)

fu_months <- as.numeric(
  sub(
    "^fu time:\\s*",
    "",
    fu_raw,
    ignore.case = TRUE
  )
)

OS_time_days <- fu_months * 30.44

geo_clinical <- data.frame(
  sample_id = rownames(pheno_geo),
  vital = vital,
  OS_event = OS_event,
  OS_time_months = fu_months,
  OS_time_days = OS_time_days,
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("GSE41613 SURVIVAL\n")
cat("========================================\n")

cat("Samples:", nrow(geo_clinical), "\n")

cat("\nVital status:\n")
print(
  table(
    geo_clinical$vital,
    useNA = "ifany"
  )
)

cat("\nOS event:\n")
print(
  table(
    geo_clinical$OS_event,
    useNA = "ifany"
  )
)

cat("\nFollow-up months:\n")
print(
  summary(
    geo_clinical$OS_time_months
  )
)

cat("\nComplete OS:\n")
print(
  sum(
    !is.na(geo_clinical$OS_event) &
    !is.na(geo_clinical$OS_time_days) &
    geo_clinical$OS_time_days > 0
  )
)

probe_ids <- rownames(expr_geo)

probe_symbol <- AnnotationDbi::mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

probe_map <- data.frame(
  probe_id = probe_ids,
  gene_name = as.character(probe_symbol),
  stringsAsFactors = FALSE
)

cat("\nMapped probes:\n")

print(
  c(
    total = nrow(probe_map),
    mapped = sum(
      !is.na(probe_map$gene_name) &
      probe_map$gene_name != ""
    )
  )
)

frozen_genes <- c(
  "CTSG",
  "P4HA1",
  "NKAIN1",
  "LRRC38",
  "SPSB4",
  "SPINK6",
  "TNFRSF4",
  "DEFB1"
)

cat("\n========================================\n")
cat("FROZEN 8-GENE PROBE COVERAGE\n")
cat("========================================\n")

probe_counts <- table(
  probe_map$gene_name
)

coverage <- data.frame(
  gene_name = frozen_genes,
  probe_n = as.integer(
    probe_counts[frozen_genes]
  ),
  stringsAsFactors = FALSE
)

coverage$probe_n[
  is.na(coverage$probe_n)
] <- 0

print(coverage)

mapped_idx <- (
  !is.na(probe_map$gene_name) &
  probe_map$gene_name != ""
)

expr_mapped <- expr_geo[
  mapped_idx,
  ,
  drop = FALSE
]

gene_mapped <- probe_map$gene_name[
  mapped_idx
]

unique_genes <- unique(
  gene_mapped
)

gene_expr_geo <- matrix(
  NA_real_,
  nrow = length(unique_genes),
  ncol = ncol(expr_mapped),
  dimnames = list(
    unique_genes,
    colnames(expr_mapped)
  )
)

for (g in unique_genes) {

  idx <- which(
    gene_mapped == g
  )

  if (length(idx) == 1) {

    gene_expr_geo[g, ] <-
      expr_mapped[idx, ]

  } else {

    gene_expr_geo[g, ] <-
      colMeans(
        expr_mapped[
          idx,
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )
  }
}

available_genes <- intersect(
  frozen_genes,
  rownames(gene_expr_geo)
)

missing_genes <- setdiff(
  frozen_genes,
  rownames(gene_expr_geo)
)

cat("\nFrozen genes available:\n")
print(available_genes)

cat("\nFrozen genes missing:\n")
print(missing_genes)

if (length(missing_genes) == 0) {

  geo_expr8 <- gene_expr_geo[
    frozen_genes,
    ,
    drop = FALSE
  ]

} else {

  geo_expr8 <- gene_expr_geo[
    available_genes,
    ,
    drop = FALSE
  ]
}

cat("\n8-gene expression dimensions:\n")
print(dim(geo_expr8))

cat("\nExpression summary:\n")
print(
  summary(
    as.numeric(geo_expr8)
  )
)

stopifnot(
  identical(
    colnames(geo_expr8),
    rownames(pheno_geo)
  )
)

saveRDS(
  geo_clinical,
  "GSE41613_OS_clinical.rds"
)

saveRDS(
  probe_map,
  "GSE41613_GPL570_probe_symbol_map.rds"
)

saveRDS(
  gene_expr_geo,
  "GSE41613_gene_expression.rds"
)

saveRDS(
  geo_expr8,
  "GSE41613_frozen8_expression.rds"
)

write.table(
  coverage,
  "GSE41613_frozen8_probe_coverage.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  geo_clinical,
  "GSE41613_OS_clinical.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Stage 7B completed.\n")
cat("========================================\n")
