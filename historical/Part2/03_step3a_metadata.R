project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

files <- c(
  "TCGA_OSCC_Nature8_raw_count_matrix.rds",
  "TCGA_OSCC_Nature8_VST_matrix.rds",
  "TCGA_OSCC_Nature8_VST_normalized_matrix.rds"
)

for(f in files){

  file_path <- file.path(proj_dir, f)

  cat("\n==============================\n")
  cat("Checking:", f, "\n")

  if(file.exists(file_path)){

    obj <- readRDS(file_path)

    cat("Class:",
        class(obj),
        "\n")

    cat("Dimension:",
        dim(obj),
        "\n")

    cat("First values:\n")

    print(
      obj[1:5,1:5]
    )

  }else{

    cat("File not found\n")

  }

}
