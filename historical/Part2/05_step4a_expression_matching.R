project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

library(DESeq2)

count_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_raw_count_matrix.rds"
)

counts_raw <- readRDS(count_file)

cat("Raw count dimension:",
    dim(counts_raw),
    "\n")

expr_sample <- colnames(counts_raw)[-(1:3)]

cat("Expression samples:",
    length(expr_sample),
    "\n")

exp_file <- file.path(
  proj_dir,
  "TCGA_OSCC_exposure_metadata_final.rds"
)

exposure <- readRDS(exp_file)

cat("Exposure rows:",
    nrow(exposure),
    "\n")

expr_patient <- substr(
  expr_sample,
  1,
  12
)

match_check <- data.frame(
  patient_id = expr_patient
)

match_check$available <-
  match_check$patient_id %in%
  exposure$patient_id

cat("\nMatching summary:\n")

print(
  table(match_check$available)
)

cat("\nSmoking table after matching:\n")

print(
  table(
    exposure$Smoking_group,
    useNA="ifany"
  )
)
