project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

setwd(project_root)

oscc_sites <- c(
  "Lip",
  "Palate",
  "Gum",
  "Base of tongue",
  "Floor of mouth",
  "Other and ill-defined sites in lip, oral cavity and pharynx",
  "Other and unspecified parts of tongue",
  "Other and unspecified parts of mouth"
)

site_patient <- unique(
  clinical[
    clinical$cases.primary_site %in% oscc_sites,
    c(
      "cases.submitter_id",
      "cases.primary_site"
    )
  ]
)

colnames(site_patient) <- c(
  "patient_barcode",
  "primary_site"
)

rna_patients <- unique(meta$patient_barcode)

site_patient <- site_patient[
  site_patient$patient_barcode %in% rna_patients,
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("Final anatomical OSCC patients with RNA-seq:\n")
cat(length(unique(site_patient$patient_barcode)), "\n")

primary_pathology <- unique(
  clinical[
    clinical$cases.submitter_id %in% site_patient$patient_barcode &
      tolower(
        as.character(
          clinical$diagnoses.diagnosis_is_primary_disease
        )
      ) == "true",
    c(
      "cases.submitter_id",
      "diagnoses.primary_diagnosis",
      "diagnoses.morphology",
      "diagnoses.icd_10_code"
    )
  ]
)

cat("\n========================================\n")
cat("PRIMARY disease diagnosis distribution:\n")

print(
  sort(
    table(
      primary_pathology$diagnoses.primary_diagnosis,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nPRIMARY disease morphology distribution:\n")

print(
  sort(
    table(
      primary_pathology$diagnoses.morphology,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

primary_non_scc <- primary_pathology[
  !grepl(
    "squamous",
    primary_pathology$diagnoses.primary_diagnosis,
    ignore.case = TRUE
  ) &
    primary_pathology$diagnoses.primary_diagnosis != "Not Reported" &
    primary_pathology$diagnoses.primary_diagnosis != "'--" &
    !is.na(primary_pathology$diagnoses.primary_diagnosis),
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("Potential non-SCC PRIMARY diagnoses:\n")

if (nrow(primary_non_scc) == 0) {

  cat("NONE\n")

} else {

  print(primary_non_scc)

}

final_meta <- meta[
  meta$patient_barcode %in% site_patient$patient_barcode,
  ,
  drop = FALSE
]

final_meta <- merge(
  final_meta,
  site_patient,
  by = "patient_barcode",
  all.x = TRUE,
  sort = FALSE
)

sample_order <- match(
  final_meta$sample_barcode,
  meta$sample_barcode
)

final_meta <- final_meta[
  order(sample_order),
  ,
  drop = FALSE
]

rownames(final_meta) <- NULL

cat("\n========================================\n")
cat("FINAL OSCC RNA cohort:\n")

cat(
  "Patients =",
  length(unique(final_meta$patient_barcode)),
  "\n"
)

cat(
  "RNA samples =",
  nrow(final_meta),
  "\n"
)

cat("\nSample types:\n")

print(
  table(
    final_meta$sample_type,
    useNA = "ifany"
  )
)

cat("\nPrimary sites × sample types:\n")

print(
  with(
    final_meta,
    table(
      primary_site,
      sample_type
    )
  )
)

cat("\n========================================\n")
cat("Duplicated sample barcodes:\n")

dup_samples <- final_meta$sample_barcode[
  duplicated(final_meta$sample_barcode)
]

if (length(dup_samples) == 0) {

  cat("NONE\n")

} else {

  print(dup_samples)

}

annotation_cols <- intersect(
  c(
    "gene_id",
    "gene_name",
    "gene_type"
  ),
  colnames(count_mat)
)

final_sample_ids <- final_meta$sample_barcode

stopifnot(
  all(final_sample_ids %in% colnames(count_mat))
)

oscc_count_all <- count_mat[
  ,
  c(
    annotation_cols,
    final_sample_ids
  ),
  drop = FALSE
]

tumor_ids <- final_meta$sample_barcode[
  final_meta$sample_type == "Primary Tumor"
]

normal_ids <- final_meta$sample_barcode[
  final_meta$sample_type == "Solid Tissue Normal"
]

oscc_count_tumor <- count_mat[
  ,
  c(
    annotation_cols,
    tumor_ids
  ),
  drop = FALSE
]

oscc_count_normal <- count_mat[
  ,
  c(
    annotation_cols,
    normal_ids
  ),
  drop = FALSE
]

cat("\n========================================\n")
cat("Final matrix dimensions:\n")

cat(
  "All OSCC:",
  nrow(oscc_count_all),
  "genes x",
  length(final_sample_ids),
  "samples\n"
)

cat(
  "Tumor:",
  nrow(oscc_count_tumor),
  "genes x",
  length(tumor_ids),
  "samples\n"
)

cat(
  "Normal:",
  nrow(oscc_count_normal),
  "genes x",
  length(normal_ids),
  "samples\n"
)

stopifnot(
  identical(
    colnames(oscc_count_all)[-(1:length(annotation_cols))],
    final_meta$sample_barcode
  )
)

stopifnot(
  length(tumor_ids) ==
    sum(final_meta$sample_type == "Primary Tumor")
)

stopifnot(
  length(normal_ids) ==
    sum(final_meta$sample_type == "Solid Tissue Normal")
)

saveRDS(
  final_meta,
  "TCGA_OSCC_Nature8_metadata.rds"
)

saveRDS(
  oscc_count_all,
  "TCGA_OSCC_Nature8_raw_count_matrix.rds"
)

saveRDS(
  oscc_count_tumor,
  "TCGA_OSCC_Nature8_tumor_raw_count_matrix.rds"
)

saveRDS(
  oscc_count_normal,
  "TCGA_OSCC_Nature8_normal_raw_count_matrix.rds"
)

saveRDS(
  primary_pathology,
  "TCGA_OSCC_Nature8_primary_pathology.rds"
)

write.table(
  final_meta,
  "TCGA_OSCC_Nature8_metadata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  primary_pathology,
  "TCGA_OSCC_Nature8_primary_pathology.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Stage 1B-3 completed successfully.\n")
cat("Nature8 OSCC cohort has been frozen.\n")
cat("========================================\n")
