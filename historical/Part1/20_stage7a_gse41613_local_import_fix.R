rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(GEOquery)

gse_file <- "GSE41613_series_matrix.txt.gz"

stopifnot(
  file.exists(gse_file)
)

gse_local <- GEOquery::getGEO(
  filename = gse_file,
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (is.list(gse_local)) {

  eset <- gse_local[[1]]

} else {

  eset <- gse_local
}

expr_geo <- Biobase::exprs(
  eset
)

pheno_geo <- Biobase::pData(
  eset
)

cat("\n========================================\n")
cat("GSE41613 LOCAL MATRIX\n")
cat("========================================\n")

cat("\nExpression dimensions:\n")
print(
  dim(expr_geo)
)

cat("\nExpression summary:\n")
print(
  summary(
    as.numeric(expr_geo)
  )
)

cat("\nPhenotype dimensions:\n")
print(
  dim(pheno_geo)
)

cat("\nPhenotype columns:\n")
print(
  colnames(pheno_geo)
)

char_cols <- grep(
  "characteristics",
  colnames(pheno_geo),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nCharacteristics columns:\n")
print(
  char_cols
)

for (x in char_cols) {

  cat("\n--------------------------------\n")
  cat(x, "\n")

  print(
    head(
      pheno_geo[[x]],
      12
    )
  )
}

surv_cols <- grep(
  "surv|death|dead|alive|follow|status|outcome|time",
  colnames(pheno_geo),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nPotential survival columns:\n")

print(
  surv_cols
)

basic_cols <- intersect(
  c(
    "title",
    "geo_accession",
    "source_name_ch1"
  ),
  colnames(pheno_geo)
)

cat("\nBasic sample information:\n")

print(
  head(
    pheno_geo[
      ,
      basic_cols,
      drop = FALSE
    ],
    10
  )
)

saveRDS(
  eset,
  "GSE41613_ExpressionSet_local.rds"
)

saveRDS(
  expr_geo,
  "GSE41613_expression_local.rds"
)

saveRDS(
  pheno_geo,
  "GSE41613_phenotype_local.rds"
)

cat("\n========================================\n")
cat("Stage 7A-local completed.\n")
cat("GPL570 was NOT downloaded.\n")
cat("========================================\n")
