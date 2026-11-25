project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

setwd(project_root)

final_meta <- readRDS(
  "TCGA_OSCC_Nature8_metadata.rds"
)

tumor_meta <- final_meta[
  final_meta$sample_type == "Primary Tumor",
  ,
  drop = FALSE
]

tumor_patients <- unique(
  tumor_meta$patient_barcode
)

cat("\n========================================\n")
cat("Frozen OSCC tumor patients:\n")
cat(length(tumor_patients), "\n")

surv_raw <- clinical[
  clinical$cases.submitter_id %in% tumor_patients,
  c(
    "cases.submitter_id",
    "demographic.vital_status",
    "demographic.days_to_death",
    "diagnoses.days_to_last_follow_up"
  ),
  drop = FALSE
]

colnames(surv_raw) <- c(
  "patient_barcode",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up"
)

clean_na <- function(x) {

  x <- as.character(x)

  x[
    x %in% c(
      "'--",
      "--",
      "",
      "Not Reported",
      "not reported",
      "Unknown",
      "unknown"
    )
  ] <- NA

  return(x)
}

surv_raw$vital_status <-
  clean_na(surv_raw$vital_status)

surv_raw$days_to_death <-
  suppressWarnings(
    as.numeric(
      clean_na(surv_raw$days_to_death)
    )
  )

surv_raw$days_to_last_follow_up <-
  suppressWarnings(
    as.numeric(
      clean_na(surv_raw$days_to_last_follow_up)
    )
  )

get_max_non_na <- function(x) {

  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  max(x)
}

get_first_non_na <- function(x) {

  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  unique(x)[1]
}

patient_ids <- unique(surv_raw$patient_barcode)

survival_patient <- data.frame(
  patient_barcode = patient_ids,
  stringsAsFactors = FALSE
)

survival_patient$vital_status <- sapply(
  patient_ids,
  function(id) {

    get_first_non_na(
      surv_raw$vital_status[
        surv_raw$patient_barcode == id
      ]
    )

  }
)

survival_patient$days_to_death <- sapply(
  patient_ids,
  function(id) {

    get_max_non_na(
      surv_raw$days_to_death[
        surv_raw$patient_barcode == id
      ]
    )

  }
)

survival_patient$days_to_last_follow_up <- sapply(
  patient_ids,
  function(id) {

    get_max_non_na(
      surv_raw$days_to_last_follow_up[
        surv_raw$patient_barcode == id
      ]
    )

  }
)

survival_patient$OS_event <- ifelse(
  tolower(survival_patient$vital_status) == "dead",
  1,
  ifelse(
    tolower(survival_patient$vital_status) == "alive",
    0,
    NA
  )
)

survival_patient$OS_time <- ifelse(
  survival_patient$OS_event == 1,
  survival_patient$days_to_death,
  survival_patient$days_to_last_follow_up
)

demo_raw <- clinical[
  clinical$cases.submitter_id %in% tumor_patients,
  c(
    "cases.submitter_id",
    "demographic.age_at_index",
    "demographic.sex_at_birth"
  ),
  drop = FALSE
]

colnames(demo_raw) <- c(
  "patient_barcode",
  "age",
  "sex"
)

demo_raw$age <- suppressWarnings(
  as.numeric(
    clean_na(demo_raw$age)
  )
)

demo_raw$sex <- clean_na(
  demo_raw$sex
)

demo_patient <- data.frame(
  patient_barcode = tumor_patients,
  stringsAsFactors = FALSE
)

demo_patient$age <- sapply(
  tumor_patients,
  function(id) {

    x <- demo_raw$age[
      demo_raw$patient_barcode == id
    ]

    x <- x[!is.na(x)]

    if (length(x) == 0) {
      return(NA_real_)
    }

    unique(x)[1]

  }
)

demo_patient$sex <- sapply(
  tumor_patients,
  function(id) {

    get_first_non_na(
      demo_raw$sex[
        demo_raw$patient_barcode == id
      ]
    )

  }
)

dx_raw <- clinical[
  clinical$cases.submitter_id %in% tumor_patients &
    tolower(
      as.character(
        clinical$diagnoses.diagnosis_is_primary_disease
      )
    ) == "true",
  c(
    "cases.submitter_id",
    "diagnoses.ajcc_pathologic_stage",
    "diagnoses.ajcc_pathologic_t",
    "diagnoses.ajcc_pathologic_n",
    "diagnoses.ajcc_pathologic_m"
  ),
  drop = FALSE
]

colnames(dx_raw) <- c(
  "patient_barcode",
  "stage",
  "T",
  "N",
  "M"
)

dx_raw$stage <- clean_na(dx_raw$stage)
dx_raw$T     <- clean_na(dx_raw$T)
dx_raw$N     <- clean_na(dx_raw$N)
dx_raw$M     <- clean_na(dx_raw$M)

count_unique_non_na <- function(x) {

  length(
    unique(
      x[!is.na(x)]
    )
  )

}

conflict_check <- data.frame(
  patient_barcode = tumor_patients,
  stage_n = NA_integer_,
  T_n = NA_integer_,
  N_n = NA_integer_,
  M_n = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(tumor_patients)) {

  id <- tumor_patients[i]

  tmp <- dx_raw[
    dx_raw$patient_barcode == id,
    ,
    drop = FALSE
  ]

  conflict_check$stage_n[i] <-
    count_unique_non_na(tmp$stage)

  conflict_check$T_n[i] <-
    count_unique_non_na(tmp$T)

  conflict_check$N_n[i] <-
    count_unique_non_na(tmp$N)

  conflict_check$M_n[i] <-
    count_unique_non_na(tmp$M)
}

conflicts <- conflict_check[
  conflict_check$stage_n > 1 |
    conflict_check$T_n > 1 |
    conflict_check$N_n > 1 |
    conflict_check$M_n > 1,
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("Patients with conflicting primary AJCC records:\n")

if (nrow(conflicts) == 0) {

  cat("NONE\n")

} else {

  print(conflicts)

}

dx_patient <- data.frame(
  patient_barcode = tumor_patients,
  stringsAsFactors = FALSE
)

for (v in c("stage", "T", "N", "M")) {

  dx_patient[[v]] <- sapply(
    tumor_patients,
    function(id) {

      get_first_non_na(
        dx_raw[
          dx_raw$patient_barcode == id,
          v
        ]
      )

    }
  )

}

clinical_model <- merge(
  survival_patient,
  demo_patient,
  by = "patient_barcode",
  all = TRUE
)

clinical_model <- merge(
  clinical_model,
  dx_patient,
  by = "patient_barcode",
  all = TRUE
)

site_info <- unique(
  tumor_meta[
    ,
    c(
      "patient_barcode",
      "primary_site"
    )
  ]
)

clinical_model <- merge(
  clinical_model,
  site_info,
  by = "patient_barcode",
  all.x = TRUE
)

clinical_model <- clinical_model[
  match(
    tumor_meta$patient_barcode,
    clinical_model$patient_barcode
  ),
  ,
  drop = FALSE
]

rownames(clinical_model) <- NULL

cat("\n========================================\n")
cat("Patient-level clinical table:\n")

cat(
  "Rows =",
  nrow(clinical_model),
  "\n"
)

cat(
  "Unique patients =",
  length(unique(clinical_model$patient_barcode)),
  "\n"
)

cat("\nVital status:\n")
print(
  table(
    clinical_model$vital_status,
    useNA = "ifany"
  )
)

cat("\nOS event:\n")
print(
  table(
    clinical_model$OS_event,
    useNA = "ifany"
  )
)

cat("\nOS time summary (days):\n")
print(
  summary(
    clinical_model$OS_time
  )
)

cat("\nNumber with complete OS information:\n")

complete_os <- (
  !is.na(clinical_model$OS_event) &
  !is.na(clinical_model$OS_time) &
  clinical_model$OS_time > 0
)

cat(
  sum(complete_os),
  "/",
  nrow(clinical_model),
  "\n"
)

cat("\n========================================\n")
cat("Age summary:\n")
print(summary(clinical_model$age))

cat("\nSex:\n")
print(
  table(
    clinical_model$sex,
    useNA = "ifany"
  )
)

cat("\nAJCC pathological stage:\n")
print(
  table(
    clinical_model$stage,
    useNA = "ifany"
  )
)

cat("\nPathological T:\n")
print(
  table(
    clinical_model$T,
    useNA = "ifany"
  )
)

cat("\nPathological N:\n")
print(
  table(
    clinical_model$N,
    useNA = "ifany"
  )
)

cat("\nPathological M:\n")
print(
  table(
    clinical_model$M,
    useNA = "ifany"
  )
)

variables_to_check <- c(
  "OS_event",
  "OS_time",
  "age",
  "sex",
  "stage",
  "T",
  "N",
  "M",
  "primary_site"
)

missing_summary <- data.frame(
  variable = variables_to_check,
  missing_n = sapply(
    clinical_model[
      variables_to_check
    ],
    function(x) sum(is.na(x))
  ),
  stringsAsFactors = FALSE
)

missing_summary$missing_percent <-
  round(
    missing_summary$missing_n /
      nrow(clinical_model) * 100,
    2
  )

cat("\n========================================\n")
cat("Missing data summary:\n")
print(missing_summary)

saveRDS(
  clinical_model,
  "TCGA_OSCC_Nature8_clinical_model_raw.rds"
)

write.table(
  clinical_model,
  "TCGA_OSCC_Nature8_clinical_model_raw.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  missing_summary,
  "TCGA_OSCC_Nature8_clinical_missing_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Stage 2A completed.\n")
cat("No patients have been removed for survival yet.\n")
cat("No Cox model has been fitted yet.\n")
cat("========================================\n")
