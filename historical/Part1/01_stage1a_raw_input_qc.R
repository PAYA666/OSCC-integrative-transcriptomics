rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

count_file <- "TCGA_HNSC_clean_count_matrix.rds"
meta_file  <- "TCGA_HNSC_sample_metadata.rds"

count_mat <- readRDS(count_file)
meta <- readRDS(meta_file)

cat("\n==============================\n")
cat("Count matrix class:\n")
print(class(count_mat))

cat("\nCount matrix dimensions:\n")
print(dim(count_mat))

cat("\nMetadata class:\n")
print(class(meta))

cat("\nMetadata dimensions:\n")
print(dim(meta))

cat("\n==============================\n")
cat("Count matrix column names (first 20):\n")
print(head(colnames(count_mat), 20))

cat("\nCount matrix first 6 rows / first 10 columns:\n")
print(count_mat[1:min(6, nrow(count_mat)),
                1:min(10, ncol(count_mat))])

cat("\n==============================\n")
cat("Metadata column names:\n")
print(colnames(meta))

cat("\nMetadata first 6 rows:\n")
print(head(meta))

possible_cols <- c(
  "sample_id",
  "sample_submitter_id",
  "case_id",
  "case_submitter_id",
  "patient_barcode",
  "sample_barcode",
  "sample_type",
  "primary_site",
  "tissue_or_organ_of_origin",
  "site_of_resection_or_biopsy",
  "primary_diagnosis",
  "morphology",
  "project_id"
)

existing_cols <- intersect(possible_cols, colnames(meta))

cat("\n==============================\n")
cat("Potentially useful metadata columns found:\n")
print(existing_cols)

for (x in existing_cols) {

  cat("\n---------------------------------\n")
  cat("Variable:", x, "\n")

  tmp <- meta[[x]]

  print(
    sort(
      table(tmp, useNA = "ifany"),
      decreasing = TRUE
    )
  )
}

annotation_cols <- intersect(
  c("gene_id", "gene_name", "gene_type"),
  colnames(count_mat)
)

expression_cols <- setdiff(
  colnames(count_mat),
  annotation_cols
)

cat("\n==============================\n")
cat("Annotation columns:\n")
print(annotation_cols)

cat("\nNumber of expression/sample columns:\n")
print(length(expression_cols))

cat("\nFirst 20 expression column names:\n")
print(head(expression_cols, 20))

cat("\n==============================\n")

candidate_id_cols <- intersect(
  c(
    "sample_id",
    "sample_submitter_id",
    "sample_barcode",
    "file_name",
    "file_id"
  ),
  colnames(meta)
)

cat("Candidate sample-ID columns in metadata:\n")
print(candidate_id_cols)

for (x in candidate_id_cols) {

  overlap_n <- sum(expression_cols %in% meta[[x]])

  cat(
    "\nDirect overlap with",
    x,
    "=",
    overlap_n,
    "/",
    length(expression_cols),
    "\n"
  )
}

cat("\n==============================\n")
cat("Stage 1A finished.\n")
cat("NO samples were filtered in this step.\n")
cat("==============================\n")
