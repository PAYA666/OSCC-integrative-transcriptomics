project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

me_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes.rds"
)

stopifnot(file.exists(me_file))

ME_old <- readRDS(me_file)

expr_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

stopifnot(file.exists(expr_file))

datExpr_old <- readRDS(expr_file)

wgcna_sample_id <- rownames(datExpr_old)
wgcna_patient_id <- substr(
  wgcna_sample_id,
  1,
  12
)

smoking_file <- file.path(
  proj_dir,
  "Part2_Smoking_sample_metadata_patient341.rds"
)

stopifnot(file.exists(smoking_file))

smoking_meta <- readRDS(smoking_file)

smoking_aligned <- smoking_meta[
  match(
    wgcna_patient_id,
    smoking_meta$patient_id
  ),
]

stopifnot(
  all(
    smoking_aligned$patient_id ==
    wgcna_patient_id
  )
)

cat("ME dimensions:", dim(ME_old), "\n")
cat("WGCNA samples:", length(wgcna_patient_id), "\n")

if (!is.null(rownames(ME_old))) {

  cat(
    "ME rownames available:",
    TRUE,
    "\n"
  )

  me_patient <- substr(
    rownames(ME_old),
    1,
    12
  )

  cat(
    "ME patient order identical:",
    identical(
      me_patient,
      wgcna_patient_id
    ),
    "\n"
  )

}

if (!is.null(rownames(ME_old))) {

  me_patient <- substr(
    rownames(ME_old),
    1,
    12
  )

  ME_use <- ME_old[
    match(
      wgcna_patient_id,
      me_patient
    ),
    ,
    drop = FALSE
  ]

} else {

  ME_use <- ME_old

}

smoking_aligned$Smoking_numeric <- NA_real_

smoking_aligned$Smoking_numeric[
  smoking_aligned$smoking_bin == "Never"
] <- 0

smoking_aligned$Smoking_numeric[
  smoking_aligned$smoking_bin == "Ever"
] <- 1

complete_idx <- !is.na(
  smoking_aligned$Smoking_numeric
)

ME_complete <- ME_use[
  complete_idx,
  ,
  drop = FALSE
]

smoking_complete <- smoking_aligned[
  complete_idx,
]

cat(
  "Complete Smoking cases:",
  nrow(smoking_complete),
  "\n"
)

print(
  table(
    smoking_complete$smoking_bin
  )
)

module_results <- lapply(
  colnames(ME_complete),
  function(module_name) {

    y <- ME_complete[, module_name]

    fit <- lm(
      y ~ smoking_complete$Smoking_numeric
    )

    coef_tab <- summary(fit)$coefficients

    beta <- coef_tab[
      "smoking_complete$Smoking_numeric",
      "Estimate"
    ]

    se <- coef_tab[
      "smoking_complete$Smoking_numeric",
      "Std. Error"
    ]

    p <- coef_tab[
      "smoking_complete$Smoking_numeric",
      "Pr(>|t|)"
    ]

    cor_val <- cor(
      y,
      smoking_complete$Smoking_numeric,
      method = "pearson"
    )

    data.frame(
      Module = module_name,
      Beta = beta,
      SE = se,
      Correlation = cor_val,
      P = p,
      stringsAsFactors = FALSE
    )
  }
)

module_results <- do.call(
  rbind,
  module_results
)

module_results$FDR <- p.adjust(
  module_results$P,
  method = "BH"
)

module_results$Direction <- ifelse(
  module_results$Beta > 0,
  "Ever_higher",
  "Never_higher"
)

module_results <- module_results[
  order(
    module_results$FDR,
    module_results$P
  ),
]

saveRDS(
  module_results,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_ModuleTrait_RECALCULATED.rds"
  )
)

write.table(
  module_results,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_ModuleTrait_RECALCULATED.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 9B finished\n")
cat("========================================\n")

print(module_results)
