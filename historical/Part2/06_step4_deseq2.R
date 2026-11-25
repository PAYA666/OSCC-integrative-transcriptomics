project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

exp_file <- file.path(
  proj_dir,
  "Part2_TCGA_OSCC_Nature8_Exposure_341_CLEAN.rds"
)

stopifnot(file.exists(exp_file))

exposure <- readRDS(exp_file)

cat("Exposure dimension:\n")
print(dim(exposure))

cat("\nColumn names:\n")
print(colnames(exposure))

cat("\nSmoking related columns:\n")

print(
grep(
"smok|tobacco|pack|cigarette",
colnames(exposure),
ignore.case = TRUE,
value = TRUE
)
)

cat("\nSmoking_bin distribution:\n")

if("smoking_bin" %in% colnames(exposure)){

  print(
    table(
      exposure$smoking_bin,
      useNA="ifany"
    )
  )

}else{

  cat("No smoking_bin column\n")

}
