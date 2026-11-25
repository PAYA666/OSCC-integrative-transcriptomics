project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

wgcna_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

stopifnot(file.exists(wgcna_file))

datExpr <- readRDS(wgcna_file)

cat("WGCNA dimensions:",
    dim(datExpr), "\n")

sample_id <- rownames(datExpr)

if (is.null(sample_id)) {
  stop("datExpr 没有rownames")
}

patient_id <- sort(
  unique(substr(sample_id,1,12))
)

cat("Tumor samples:",
    nrow(datExpr), "\n")

cat("Patient IDs:",
    length(patient_id), "\n")

bad_id <- patient_id[
  !grepl(
    "^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}$",
    patient_id
  )
]

if(length(bad_id)>0){
  print(bad_id)
}else{
  cat("TCGA barcode check passed\n")
}

write.table(
  patient_id,
  file=file.path(
    proj_dir,
    "TCGA_OSCC_Nature8_patient_IDs_for_exposure.txt"
  ),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)
