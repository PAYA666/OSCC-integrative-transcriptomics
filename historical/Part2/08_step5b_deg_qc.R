project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

count_matrix <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_count_matrix_341.rds"
  )
)

sample_info <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_sample_metadata_341.rds"
  )
)

sample_info$sample_portion <-
  substr(
    sample_info$sample_barcode,
    14,
    16
  )

sample_info$priority <-
  ifelse(
    grepl(
      "-01A$",
      sample_info$sample_barcode
    ),
    1,
    2
  )

sample_info <- sample_info[
  order(
    sample_info$patient_id,
    sample_info$priority
  ),
]

sample_info_unique <-
  sample_info[
    !duplicated(
      sample_info$patient_id
    ),
]

keep_samples <-
  sample_info_unique$sample_barcode

count_matrix_unique <-
  count_matrix[
    ,
    keep_samples
  ]

rownames(sample_info_unique) <-
  sample_info_unique$sample_barcode

saveRDS(
  count_matrix_unique,
  file.path(
    proj_dir,
    "Part2_Smoking_count_matrix_patient341.rds"
  )
)

saveRDS(
  sample_info_unique,
  file.path(
    proj_dir,
    "Part2_Smoking_sample_metadata_patient341.rds"
  )
)

cat("Finished patient-level matrix\n")
