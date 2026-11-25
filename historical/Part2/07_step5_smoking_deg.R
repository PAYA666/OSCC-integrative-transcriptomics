project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

count_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_raw_count_matrix.rds"
)

counts_raw <- readRDS(count_file)

exp_file <- file.path(
  proj_dir,
  "Part2_TCGA_OSCC_Nature8_Exposure_341_CLEAN.rds"
)

exposure <- readRDS(exp_file)

sample_barcode <- colnames(counts_raw)[-(1:3)]

patient_from_rna <- substr(
  sample_barcode,
  1,
  12
)

keep_patient <- exposure$patient_id

keep_sample_index <- patient_from_rna %in%
  keep_patient

cat(
"RNA samples matched:",
sum(keep_sample_index),
"\n"
)

count_matrix <- counts_raw[
  ,
  c(
    FALSE,
    FALSE,
    FALSE,
    keep_sample_index
  )
]

rownames(count_matrix) <-
  counts_raw$gene_id

colnames(count_matrix) <-
  sample_barcode[keep_sample_index]

sample_info <- data.frame(

  sample_barcode =
    colnames(count_matrix),

  patient_id =
    substr(
      colnames(count_matrix),
      1,
      12
    )

)

sample_info <- merge(
  sample_info,
  exposure[,c(
    "patient_id",
    "smoking_bin"
  )],
  by="patient_id",
  all.x=TRUE
)

rownames(sample_info) <-
  sample_info$sample_barcode

saveRDS(
  count_matrix,
  file.path(
    proj_dir,
    "Part2_Smoking_count_matrix_341.rds"
  )
)

saveRDS(
  sample_info,
  file.path(
    proj_dir,
    "Part2_Smoking_sample_metadata_341.rds"
  )
)

cat("Finished\n")
