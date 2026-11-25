project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

wgcna_expr_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

wgcna_network_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_network_power8.rds"
)

wgcna_me_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes.rds"
)

wgcna_gene_module_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv"
)

smoking_meta_file <- file.path(
  proj_dir,
  "Part2_Smoking_sample_metadata_patient341.rds"
)

file_check <- data.frame(
  file = c(
    "WGCNA expression",
    "WGCNA network",
    "Module eigengenes",
    "Gene-module table",
    "Smoking metadata"
  ),
  exists = c(
    file.exists(wgcna_expr_file),
    file.exists(wgcna_network_file),
    file.exists(wgcna_me_file),
    file.exists(wgcna_gene_module_file),
    file.exists(smoking_meta_file)
  )
)

cat("\n===== File check =====\n")
print(file_check)

if (!file.exists(wgcna_expr_file)) {
  stop("Cannot find old WGCNA datExpr file.")
}

if (!file.exists(smoking_meta_file)) {
  stop("Cannot find current Smoking metadata file.")
}

datExpr_old <- readRDS(
  wgcna_expr_file
)

cat("\n===== Old WGCNA datExpr =====\n")

cat(
  "Class:",
  paste(class(datExpr_old), collapse = ", "),
  "\n"
)

cat(
  "Dimensions:",
  nrow(datExpr_old),
  "samples x",
  ncol(datExpr_old),
  "genes\n"
)

old_sample_id <- rownames(datExpr_old)

if (is.null(old_sample_id)) {
  stop("Old WGCNA datExpr has no rownames.")
}

old_patient_id <- substr(
  old_sample_id,
  1,
  12
)

smoking_meta <- readRDS(
  smoking_meta_file
)

required_cols <- c(
  "sample_barcode",
  "patient_id",
  "smoking_bin"
)

missing_cols <- setdiff(
  required_cols,
  colnames(smoking_meta)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Smoking metadata missing columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

cat("\n===== Patient uniqueness =====\n")

cat(
  "Old WGCNA rows:",
  length(old_patient_id),
  "\n"
)

cat(
  "Unique old WGCNA patients:",
  length(unique(old_patient_id)),
  "\n"
)

cat(
  "Current Smoking metadata rows:",
  nrow(smoking_meta),
  "\n"
)

cat(
  "Unique current patients:",
  length(unique(smoking_meta$patient_id)),
  "\n"
)

cat(
  "Duplicated old WGCNA patients:",
  sum(duplicated(old_patient_id)),
  "\n"
)

cat(
  "Duplicated current patients:",
  sum(duplicated(smoking_meta$patient_id)),
  "\n"
)

common_patient <- intersect(
  unique(old_patient_id),
  unique(smoking_meta$patient_id)
)

old_only <- setdiff(
  unique(old_patient_id),
  unique(smoking_meta$patient_id)
)

current_only <- setdiff(
  unique(smoking_meta$patient_id),
  unique(old_patient_id)
)

cat("\n===== Cohort overlap =====\n")

cat(
  "Common patients:",
  length(common_patient),
  "\n"
)

cat(
  "Old WGCNA only:",
  length(old_only),
  "\n"
)

cat(
  "Current cohort only:",
  length(current_only),
  "\n"
)

same_patient_set <- setequal(
  unique(old_patient_id),
  unique(smoking_meta$patient_id)
)

cat(
  "Exact same patient set:",
  same_patient_set,
  "\n"
)

smoking_match <- smoking_meta[
  match(
    old_patient_id,
    smoking_meta$patient_id
  ),
]

cat("\n===== Smoking distribution in WGCNA cohort =====\n")

print(
  table(
    smoking_match$smoking_bin,
    useNA = "ifany"
  )
)

cat("\n===== Existing WGCNA objects =====\n")

if (file.exists(wgcna_network_file)) {

  network_old <- readRDS(
    wgcna_network_file
  )

  cat(
    "Network object class:",
    paste(class(network_old), collapse = ", "),
    "\n"
  )

  if (is.list(network_old)) {

    cat(
      "Network object names:\n"
    )

    print(
      names(network_old)
    )
  }
}

if (file.exists(wgcna_me_file)) {

  ME_old <- readRDS(
    wgcna_me_file
  )

  cat(
    "\nModule eigengene dimensions:",
    paste(dim(ME_old), collapse = " x "),
    "\n"
  )

  cat(
    "Number of module eigengenes:",
    ncol(ME_old),
    "\n"
  )
}

if (file.exists(wgcna_gene_module_file)) {

  gene_module_old <- read.delim(
    wgcna_gene_module_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cat(
    "\nGene-module table dimensions:",
    paste(dim(gene_module_old), collapse = " x "),
    "\n"
  )

  cat(
    "Gene-module columns:\n"
  )

  print(
    colnames(gene_module_old)
  )
}

audit_summary <- data.frame(
  Item = c(
    "Old_WGCNA_samples",
    "Old_WGCNA_unique_patients",
    "Current_unique_patients",
    "Common_patients",
    "Old_only",
    "Current_only",
    "Exact_same_patient_set",
    "Smoking_available_in_old_WGCNA"
  ),

  Value = c(
    nrow(datExpr_old),
    length(unique(old_patient_id)),
    length(unique(smoking_meta$patient_id)),
    length(common_patient),
    length(old_only),
    length(current_only),
    as.character(same_patient_set),
    sum(!is.na(smoking_match$smoking_bin))
  ),

  stringsAsFactors = FALSE
)

write.table(
  audit_summary,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_reuse_audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 9A audit finished\n")
cat("========================================\n")
