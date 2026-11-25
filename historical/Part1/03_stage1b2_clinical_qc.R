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

cat("\n========================================\n")
cat("Pre-specified OSCC anatomical sites:\n")
print(oscc_sites)

site_patient <- unique(
  clinical[, c(
    "cases.submitter_id",
    "cases.primary_site"
  )]
)

cat("\n========================================\n")
cat("Primary-site distribution at patient level:\n")

print(
  sort(
    table(
      site_patient$cases.primary_site,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

oscc_site_patient <- site_patient[
  site_patient$cases.primary_site %in% oscc_sites,
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("OSCC patients selected by anatomical site:\n")
cat(
  length(unique(oscc_site_patient$cases.submitter_id)),
  "\n"
)

cat("\nPatient number by selected primary site:\n")

print(
  sort(
    table(oscc_site_patient$cases.primary_site),
    decreasing = TRUE
  )
)

site_per_patient <- table(
  oscc_site_patient$cases.submitter_id
)

cat("\n========================================\n")
cat("Patients with >1 selected primary_site:\n")

print(
  site_per_patient[site_per_patient > 1]
)

rna_patients <- unique(meta$patient_barcode)

oscc_rna_patients <- intersect(
  unique(oscc_site_patient$cases.submitter_id),
  rna_patients
)

cat("\n========================================\n")
cat("Selected OSCC patients with RNA-seq:\n")
cat(length(oscc_rna_patients), "\n")

oscc_meta_tmp <- meta[
  meta$patient_barcode %in% oscc_rna_patients,
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("RNA-seq samples belonging to selected OSCC patients:\n")
cat(nrow(oscc_meta_tmp), "\n")

cat("\nSample type distribution:\n")

print(
  table(
    oscc_meta_tmp$sample_type,
    useNA = "ifany"
  )
)

site_map <- unique(
  oscc_site_patient[
    oscc_site_patient$cases.submitter_id %in%
      oscc_rna_patients,
    c(
      "cases.submitter_id",
      "cases.primary_site"
    )
  ]
)

colnames(site_map) <- c(
  "patient_barcode",
  "primary_site"
)

oscc_meta_check <- merge(
  oscc_meta_tmp,
  site_map,
  by = "patient_barcode",
  all.x = TRUE
)

cat("\n========================================\n")
cat("RNA sample distribution by primary site and sample type:\n")

print(
  with(
    oscc_meta_check,
    table(
      primary_site,
      sample_type,
      useNA = "ifany"
    )
  )
)

pathology <- unique(
  clinical[
    clinical$cases.submitter_id %in% oscc_rna_patients,
    c(
      "cases.submitter_id",
      "diagnoses.primary_diagnosis",
      "diagnoses.morphology",
      "diagnoses.icd_10_code",
      "diagnoses.diagnosis_is_primary_disease"
    )
  ]
)

cat("\n========================================\n")
cat("Primary diagnosis distribution:\n")

print(
  sort(
    table(
      pathology$diagnoses.primary_diagnosis,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\nMorphology distribution:\n")

print(
  sort(
    table(
      pathology$diagnoses.morphology,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

non_scc_records <- pathology[
  !grepl(
    "squamous",
    pathology$diagnoses.primary_diagnosis,
    ignore.case = TRUE
  ) &
    !is.na(pathology$diagnoses.primary_diagnosis) &
    pathology$diagnoses.primary_diagnosis != "'--",
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("Records without 'squamous' in primary diagnosis:\n")

if (nrow(non_scc_records) == 0) {

  cat("NONE\n")

} else {

  print(non_scc_records)

}

annotation_cols <- c(
  "gene_id",
  "gene_name",
  "gene_type"
)

expression_cols <- setdiff(
  colnames(count_mat),
  annotation_cols
)

matched_samples <- oscc_meta_tmp$sample_barcode %in%
  expression_cols

cat("\n========================================\n")
cat("OSCC RNA samples matched to count matrix:\n")

cat(
  sum(matched_samples),
  "/",
  nrow(oscc_meta_tmp),
  "\n"
)

if (any(!matched_samples)) {

  cat("\nUnmatched samples:\n")

  print(
    oscc_meta_tmp[
      !matched_samples,
      c(
        "patient_barcode",
        "sample_barcode",
        "sample_type"
      ),
      drop = FALSE
    ]
  )

}

write.table(
  oscc_meta_check,
  file = "Stage1B2_OSCC_candidate_metadata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pathology,
  file = "Stage1B2_OSCC_pathology_QC.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Stage 1B-2 completed.\n")
cat("No final expression matrix has been created yet.\n")
cat("Candidate metadata and pathology QC tables saved.\n")
cat("========================================\n")
