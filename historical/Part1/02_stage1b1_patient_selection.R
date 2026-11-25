project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

setwd(project_root)

clinical_file <- "clinical.tsv"

clinical <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\n==============================\n")
cat("Clinical dimensions:\n")
print(dim(clinical))

cat("\n==============================\n")
cat("All clinical column names:\n")
print(colnames(clinical))

keywords <- c(
  "submitter",
  "site",
  "tissue",
  "organ",
  "diagnosis",
  "morphology",
  "icd",
  "anatomic"
)

for (kw in keywords) {

  cat("\n---------------------------------\n")
  cat("Keyword:", kw, "\n")

  print(
    grep(
      kw,
      colnames(clinical),
      value = TRUE,
      ignore.case = TRUE
    )
  )
}

candidate_cols <- unique(
  unlist(
    lapply(
      c(
        "submitter",
        "site",
        "tissue",
        "organ",
        "diagnosis",
        "morphology",
        "icd"
      ),
      function(x) {
        grep(
          x,
          colnames(clinical),
          value = TRUE,
          ignore.case = TRUE
        )
      }
    )
  )
)

cat("\n==============================\n")
cat("Potentially relevant clinical fields:\n")
print(candidate_cols)

cat("\nFirst 10 rows of these fields:\n")

if (length(candidate_cols) > 0) {

  print(
    clinical[
      1:min(10, nrow(clinical)),
      candidate_cols,
      drop = FALSE
    ]
  )

}

cat("\n==============================\n")

if ("cases.submitter_id" %in% colnames(clinical)) {

  cat("cases.submitter_id exists.\n")

  cat(
    "Unique patients:",
    length(unique(clinical$cases.submitter_id)),
    "\n"
  )

  cat(
    "Overlap with RNA metadata:",
    sum(unique(meta$patient_barcode) %in%
          unique(clinical$cases.submitter_id)),
    "/",
    length(unique(meta$patient_barcode)),
    "\n"
  )

} else {

  cat("cases.submitter_id NOT found.\n")

}

cat("\nStage 1B-1 finished. No samples filtered.\n")
