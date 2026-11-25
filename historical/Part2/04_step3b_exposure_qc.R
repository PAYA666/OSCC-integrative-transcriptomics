project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

meta_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_metadata.tsv"
)

stopifnot(file.exists(meta_file))

meta <- read.delim(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Metadata dimension:",
    dim(meta),
    "\n\n")

cat("Column names:\n")

print(
  colnames(meta)
)

cat("\nFirst rows:\n")

print(
  meta[1:5,1:min(10,ncol(meta))]
)
