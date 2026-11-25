project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

input_file <- file.path(
  proj_dir,
  "nature8_exposure.rds"
)

stopifnot(file.exists(input_file))

nature8_exposure <- readRDS(input_file)

cat("Original dimensions:",
    dim(nature8_exposure),
    "\n")

nature8_exposure$Smoking_group <- NA

nature8_exposure$Smoking_group[
  nature8_exposure$tobacco_smoking_status ==
    "Lifelong Non-Smoker"
] <- "Non-smoker"

nature8_exposure$Smoking_group[
  grepl(
    "Current|Reformed",
    nature8_exposure$tobacco_smoking_status,
    ignore.case = TRUE
  )
] <- "Smoker"

nature8_exposure$Alcohol_group <- NA

nature8_exposure$Alcohol_group[
  nature8_exposure$alcohol_history == "Yes"
] <- "Alcohol"

nature8_exposure$Alcohol_group[
  nature8_exposure$alcohol_history == "No"
] <- "No_Alcohol"

final_exposure <- nature8_exposure[
  ,
  c(
    "patient_id",
    "Smoking_group",
    "Alcohol_group"
  )
]

saveRDS(
  final_exposure,
  file.path(
    proj_dir,
    "TCGA_OSCC_exposure_metadata_final.rds"
  )
)

write.table(
  final_exposure,
  file.path(
    proj_dir,
    "TCGA_OSCC_exposure_metadata_final.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("Saved final exposure metadata\n")
