rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required.")
}

suppressPackageStartupMessages({
  library(readxl)
})

gdsc_dir <- file.path(project_root, "Part3_Cisplatin/GDSC")

model_file <- file.path(
  gdsc_dir,
  "model_list_20260811.csv"
)

drug_file <- file.path(
  gdsc_dir,
  "GDSC1_fitted_dose_response_27Oct23.xlsx"
)

compound_file <- file.path(
  gdsc_dir,
  "screened_compounds_rel_8.5.csv"
)

tpm_file <- file.path(
  gdsc_dir,
  "rnaseq_merged_rsem_tpm_20260323.csv"
)

step6_file <- file.path(
  file.path(project_root, "Part3_Cisplatin"),
  "Step6_Candidate_Consolidation",
  "Part3_Step6_CandidateConsolidation_FROZEN.rds"
)

out_dir <- file.path(
  gdsc_dir,
  "Step7B_GDSC1_Pharmacogenomic_Validation"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  model_file,
  drug_file,
  compound_file,
  tpm_file,
  step6_file
)

if (!all(file.exists(required_files))) {

  cat("\nMissing files:\n")

  print(
    required_files[
      !file.exists(required_files)
    ]
  )

  stop("Required input file missing.")
}

ultra7 <- c(
  "CCN1",
  "JUN",
  "CCL2",
  "ZBTB10",
  "CCL20",
  "DAB2",
  "KYNU"
)

primary4 <- c(
  "CCN1",
  "CCL2",
  "ZBTB10",
  "DAB2"
)

step6 <- readRDS(
  step6_file
)

if (!"full_priority_evidence" %in% names(step6)) {
  stop("Step6 RDS does not contain 'full_priority_evidence'.")
}

evidence <- as.data.frame(
  step6[["full_priority_evidence"]]
)

required_evidence_cols <- c(
  "Gene",
  "Meta_Z",
  "In_Priority52"
)

if (!all(required_evidence_cols %in% colnames(evidence))) {

  stop(
    paste(
      "Step6 evidence missing:",
      paste(
        setdiff(
          required_evidence_cols,
          colnames(evidence)
        ),
        collapse = ", "
      )
    )
  )
}

evidence$Gene <- as.character(
  evidence$Gene
)

evidence$Meta_Z <- as.numeric(
  evidence$Meta_Z
)

priority52 <- unique(
  evidence$Gene[
    evidence$In_Priority52 %in% TRUE
  ]
)

direction_table <- evidence[
  !duplicated(evidence$Gene),
  c(
    "Gene",
    "Meta_Z"
  ),
  drop = FALSE
]

if (!all(ultra7 %in% direction_table$Gene)) {
  stop("One or more Ultra7 genes are absent from frozen Step6 evidence.")
}

extract_sidm <- function(x) {

  x <- toupper(
    trimws(
      as.character(x)
    )
  )

  hit <- regexpr(
    "SIDM[0-9]+",
    x,
    perl = TRUE
  )

  out <- rep(
    NA_character_,
    length(x)
  )

  good <- hit > 0

  if (any(good)) {
    out[good] <- regmatches(x, hit)[good]
  }

  return(out)
}

model <- read.csv(
  model_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

model$SIDM <- extract_sidm(
  model$model_id
)

hnscc <- model[
  trimws(model$tissue) == "Head and Neck" &
    trimws(model$cancer_type) %in%
    c(
      "Head and Neck Carcinoma",
      "Oral Cavity Carcinoma"
    ),
  ,
  drop = FALSE
]

hnscc <- hnscc[
  !is.na(hnscc$SIDM) &
    !duplicated(hnscc$SIDM),
  ,
  drop = FALSE
]

compound <- read.csv(
  compound_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cis_annotation <- compound[
  compound$DRUG_ID == 1005,
  ,
  drop = FALSE
]

if (nrow(cis_annotation) == 0) {
  stop("DRUG_ID 1005 absent from compound annotation.")
}

if (!any(
  grepl(
    "cisplatin",
    cis_annotation$DRUG_NAME,
    ignore.case = TRUE
  )
)) {
  stop("DRUG_ID 1005 is not annotated as Cisplatin.")
}

drug <- readxl::read_excel(
  drug_file
)

required_drug_cols <- c(
  "SANGER_MODEL_ID",
  "CELL_LINE_NAME",
  "DRUG_ID",
  "DRUG_NAME",
  "LN_IC50",
  "AUC"
)

if (!all(required_drug_cols %in% colnames(drug))) {
  stop("Unexpected GDSC1 fitted-dose-response structure.")
}

cis <- as.data.frame(
  drug[
    drug$DRUG_ID == 1005 &
      !is.na(drug$LN_IC50),
    required_drug_cols
  ]
)

cis$SIDM <- extract_sidm(
  cis$SANGER_MODEL_ID
)

cis$LN_IC50 <- as.numeric(
  cis$LN_IC50
)

cis$AUC <- as.numeric(
  cis$AUC
)

cis <- cis[
  !is.na(cis$SIDM) &
    is.finite(cis$LN_IC50),
  ,
  drop = FALSE
]

duplicated_cis <- unique(
  cis$SIDM[
    duplicated(cis$SIDM)
  ]
)

if (length(duplicated_cis) > 0) {

  cat("\nDuplicated Cisplatin models:\n")

  print(duplicated_cis)

  stop(
    "Duplicate DRUG_ID 1005 fitted-response entries found."
  )
}

cat("\nReading RNA-seq TPM with base R...\n")

tpm_raw <- read.csv(
  tpm_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "TPM raw dimension:",
  nrow(tpm_raw),
  "x",
  ncol(tpm_raw),
  "\n"
)

if (
  nrow(tpm_raw) < 100 ||
  ncol(tpm_raw) < 100
) {
  stop("Unexpected TPM matrix dimensions.")
}

rna_raw_ids <- colnames(
  tpm_raw
)[
  4:ncol(tpm_raw)
]

rna_ids <- extract_sidm(
  rna_raw_ids
)

if (any(is.na(rna_ids))) {

  bad_headers <- rna_raw_ids[
    is.na(rna_ids)
  ]

  cat("\nUnrecognized RNA model headers:\n")

  print(
    head(
      bad_headers,
      20
    )
  )

  stop("Some RNA expression columns do not contain SIDM IDs.")
}

if (anyDuplicated(rna_ids) > 0) {
  stop("Duplicated SIDM IDs in RNA expression columns.")
}

gene_rows <- 4:nrow(
  tpm_raw
)

gene_symbol <- trimws(
  as.character(
    tpm_raw[
      gene_rows,
      1
    ]
  )
)

all_targets <- unique(
  c(
    priority52,
    ultra7,
    primary4
  )
)

target_position <- which(
  gene_symbol %in% all_targets
)

if (length(target_position) == 0) {
  stop("No frozen candidate genes found in TPM matrix.")
}

target_data_rows <- gene_rows[
  target_position
]

target_symbols <- gene_symbol[
  target_position
]

target_expr <- as.matrix(
  tpm_raw[
    target_data_rows,
    4:ncol(tpm_raw),
    drop = FALSE
  ]
)

suppressWarnings(
  storage.mode(target_expr) <- "numeric"
)

rownames(target_expr) <- target_symbols

colnames(target_expr) <- rna_ids

rm(tpm_raw)

gc()

duplicate_symbols <- unique(
  rownames(target_expr)[
    duplicated(rownames(target_expr)) |
      duplicated(
        rownames(target_expr),
        fromLast = TRUE
      )
  ]
)

if (length(duplicate_symbols) > 0) {

  cat(
    "\nDuplicate frozen gene symbols collapsed by mean TPM:\n"
  )

  print(duplicate_symbols)
}

expr_sum <- rowsum(
  target_expr,
  group = rownames(target_expr),
  reorder = FALSE
)

expr_n <- table(
  rownames(target_expr)
)

target_expr <- sweep(
  expr_sum,
  1,
  as.numeric(
    expr_n[
      rownames(expr_sum)
    ]
  ),
  FUN = "/"
)

ultra_present <- intersect(
  ultra7,
  rownames(target_expr)
)

primary_present <- intersect(
  primary4,
  rownames(target_expr)
)

priority_present <- intersect(
  priority52,
  rownames(target_expr)
)

ultra_missing <- setdiff(
  ultra7,
  ultra_present
)

primary_missing <- setdiff(
  primary4,
  primary_present
)

priority_missing <- setdiff(
  priority52,
  priority_present
)

if (length(ultra_missing) > 0) {

  stop(
    paste(
      "Missing Ultra7:",
      paste(
        ultra_missing,
        collapse = ", "
      )
    )
  )
}

if (length(primary_missing) > 0) {

  stop(
    paste(
      "Missing Primary4:",
      paste(
        primary_missing,
        collapse = ", "
      )
    )
  )
}

hnscc_ids <- unique(
  hnscc$SIDM
)

cis_ids <- unique(
  cis$SIDM
)

rna_ids <- unique(
  colnames(target_expr)
)

hnscc_cis_ids <- intersect(
  hnscc_ids,
  cis_ids
)

hnscc_rna_ids <- intersect(
  hnscc_ids,
  rna_ids
)

cis_rna_ids <- intersect(
  cis_ids,
  rna_ids
)

matched_ids <- Reduce(
  intersect,
  list(
    hnscc_ids,
    cis_ids,
    rna_ids
  )
)

if (length(matched_ids) < 10) {

  stop(
    paste(
      "Unexpected triple overlap:",
      length(matched_ids)
    )
  )
}

expr_tpm <- target_expr[
  ,
  matched_ids,
  drop = FALSE
]

expr_log <- log2(
  expr_tpm + 1
)

gene_sd <- apply(
  expr_log,
  1,
  sd,
  na.rm = TRUE
)

valid_genes <- names(
  gene_sd
)[
  is.finite(gene_sd) &
    gene_sd > 0
]

if (!all(
  ultra7 %in% valid_genes
)) {

  stop(
    paste(
      "Zero-variance Ultra7 gene(s):",
      paste(
        setdiff(
          ultra7,
          valid_genes
        ),
        collapse = ", "
      )
    )
  )
}

if (!all(
  primary4 %in% valid_genes
)) {

  stop(
    paste(
      "Zero-variance Primary4 gene(s):",
      paste(
        setdiff(
          primary4,
          valid_genes
        ),
        collapse = ", "
      )
    )
  )
}

priority_valid <- intersect(
  priority_present,
  valid_genes
)

make_frozen_score <- function(
    expression_matrix,
    genes,
    direction_table
) {

  genes <- intersect(
    genes,
    rownames(expression_matrix)
  )

  x <- expression_matrix[
    genes,
    ,
    drop = FALSE
  ]

  z <- t(
    scale(
      t(x)
    )
  )

  meta_z <- direction_table$Meta_Z[
    match(
      rownames(z),
      direction_table$Gene
    )
  ]

  if (any(is.na(meta_z))) {

    missing_direction <- rownames(z)[
      is.na(meta_z)
    ]

    stop(
      paste(
        "Missing frozen Meta_Z for:",
        paste(
          missing_direction,
          collapse = ", "
        )
      )
    )
  }

  direction_sign <- ifelse(
    meta_z > 0,
    1,
    -1
  )

  signed_z <- sweep(
    z,
    1,
    direction_sign,
    FUN = "*"
  )

  score <- colMeans(
    signed_z,
    na.rm = TRUE
  )

  names(score) <- colnames(
    signed_z
  )

  return(
    list(
      score = score,
      genes_used = rownames(signed_z),
      signed_z = signed_z
    )
  )
}

ultra_score_object <- make_frozen_score(
  expression_matrix = expr_log,
  genes = ultra7,
  direction_table = direction_table
)

primary_score_object <- make_frozen_score(
  expression_matrix = expr_log,
  genes = primary4,
  direction_table = direction_table
)

priority_score_object <- make_frozen_score(
  expression_matrix = expr_log,
  genes = priority_valid,
  direction_table = direction_table
)

final <- data.frame(
  SIDM = matched_ids,
  stringsAsFactors = FALSE
)

idx_model <- match(
  final$SIDM,
  hnscc$SIDM
)

idx_cis <- match(
  final$SIDM,
  cis$SIDM
)

final$SANGER_MODEL_ID <- hnscc$model_id[
  idx_model
]

final$model_name <- hnscc$model_name[
  idx_model
]

final$cancer_type <- hnscc$cancer_type[
  idx_model
]

final$sample_site <- hnscc$sample_site[
  idx_model
]

final$LN_IC50 <- cis$LN_IC50[
  idx_cis
]

final$AUC <- cis$AUC[
  idx_cis
]

final$Ultra7_Score <- as.numeric(
  ultra_score_object$score[
    final$SIDM
  ]
)

final$Primary4_Score <- as.numeric(
  primary_score_object$score[
    final$SIDM
  ]
)

final$Priority52_Score <- as.numeric(
  priority_score_object$score[
    final$SIDM
  ]
)

for (gene in ultra7) {

  final[, gene] <- as.numeric(
    expr_log[
      gene,
      final$SIDM
    ]
  )
}

check_numeric <- c(
  "LN_IC50",
  "Ultra7_Score",
  "Primary4_Score",
  "Priority52_Score",
  ultra7
)

for (nm in check_numeric) {

  if (any(
    !is.finite(
      final[, nm]
    )
  )) {

    stop(
      paste(
        "Non-finite value detected in:",
        nm
      )
    )
  }
}

primary_test <- suppressWarnings(
  cor.test(
    final$Ultra7_Score,
    final$LN_IC50,
    method = "spearman",
    exact = FALSE
  )
)

primary_result <- data.frame(
  Analysis = "Primary",
  Cohort = "HNSCC",
  Predictor = "Ultra7_Score",
  N = nrow(final),
  Spearman_Rho = unname(
    primary_test$estimate
  ),
  P_value = primary_test$p.value,
  Expected_Direction = "Positive",
  Direction_Supported = unname(
    primary_test$estimate
  ) > 0,
  stringsAsFactors = FALSE
)

secondary_primary4 <- data.frame()

for (gene in primary4) {

  ct <- suppressWarnings(
    cor.test(
      final[, gene],
      final$LN_IC50,
      method = "spearman",
      exact = FALSE
    )
  )

  gene_meta_z <- direction_table$Meta_Z[
    match(
      gene,
      direction_table$Gene
    )
  ]

  expected_direction <- ifelse(
    gene_meta_z > 0,
    "Positive",
    "Negative"
  )

  direction_supported <- ifelse(
    gene_meta_z > 0,
    unname(ct$estimate) > 0,
    unname(ct$estimate) < 0
  )

  secondary_primary4 <- rbind(
    secondary_primary4,
    data.frame(
      Gene = gene,
      N = nrow(final),
      Spearman_Rho = unname(
        ct$estimate
      ),
      P_value = ct$p.value,
      Expected_Direction = expected_direction,
      Direction_Supported = direction_supported,
      stringsAsFactors = FALSE
    )
  )
}

secondary_primary4$FDR_BH4 <- p.adjust(
  secondary_primary4$P_value,
  method = "BH"
)

secondary_primary4 <- secondary_primary4[
  match(
    primary4,
    secondary_primary4$Gene
  ),
  ,
  drop = FALSE
]

ultra7_single <- data.frame()

for (gene in ultra7) {

  ct <- suppressWarnings(
    cor.test(
      final[, gene],
      final$LN_IC50,
      method = "spearman",
      exact = FALSE
    )
  )

  gene_meta_z <- direction_table$Meta_Z[
    match(
      gene,
      direction_table$Gene
    )
  ]

  ultra7_single <- rbind(
    ultra7_single,
    data.frame(
      Gene = gene,
      N = nrow(final),
      Spearman_Rho = unname(
        ct$estimate
      ),
      P_value = ct$p.value,
      Expected_Direction = ifelse(
        gene_meta_z > 0,
        "Positive",
        "Negative"
      ),
      Direction_Supported = ifelse(
        gene_meta_z > 0,
        unname(ct$estimate) > 0,
        unname(ct$estimate) < 0
      ),
      stringsAsFactors = FALSE
    )
  )
}

ultra7_single$FDR_BH7 <- p.adjust(
  ultra7_single$P_value,
  method = "BH"
)

ultra7_single <- ultra7_single[
  match(
    ultra7,
    ultra7_single$Gene
  ),
  ,
  drop = FALSE
]

score_sensitivity <- data.frame()

for (
  predictor in c(
    "Primary4_Score",
    "Priority52_Score"
  )
) {

  ct <- suppressWarnings(
    cor.test(
      final[, predictor],
      final$LN_IC50,
      method = "spearman",
      exact = FALSE
    )
  )

  score_sensitivity <- rbind(
    score_sensitivity,
    data.frame(
      Cohort = "HNSCC",
      Predictor = predictor,
      N = nrow(final),
      Spearman_Rho = unname(
        ct$estimate
      ),
      P_value = ct$p.value,
      Expected_Direction = "Positive",
      Direction_Supported = unname(
        ct$estimate
      ) > 0,
      stringsAsFactors = FALSE
    )
  )
}

oral_final <- final[
  trimws(final$cancer_type) ==
    "Oral Cavity Carcinoma",
  ,
  drop = FALSE
]

oral_results <- data.frame()

if (nrow(oral_final) >= 5) {

  for (
    predictor in c(
      "Ultra7_Score",
      "Primary4_Score",
      "Priority52_Score"
    )
  ) {

    ct <- suppressWarnings(
      cor.test(
        oral_final[, predictor],
        oral_final$LN_IC50,
        method = "spearman",
        exact = FALSE
      )
    )

    oral_results <- rbind(
      oral_results,
      data.frame(
        Cohort = "Oral Cavity Carcinoma",
        Predictor = predictor,
        N = nrow(oral_final),
        Spearman_Rho = unname(
          ct$estimate
        ),
        P_value = ct$p.value,
        Expected_Direction = "Positive",
        Direction_Supported = unname(
          ct$estimate
        ) > 0,
        stringsAsFactors = FALSE
      )
    )
  }

} else {

  oral_results <- data.frame(
    Cohort = "Oral Cavity Carcinoma",
    Predictor = NA_character_,
    N = nrow(oral_final),
    Spearman_Rho = NA_real_,
    P_value = NA_real_,
    Expected_Direction = NA_character_,
    Direction_Supported = NA,
    stringsAsFactors = FALSE
  )
}

qc_summary <- data.frame(
  HNSCC_Annotated = length(
    hnscc_ids
  ),

  Cisplatin_All = length(
    cis_ids
  ),

  RNA_All = length(
    rna_ids
  ),

  HNSCC_Cisplatin = length(
    hnscc_cis_ids
  ),

  HNSCC_RNA = length(
    hnscc_rna_ids
  ),

  Cisplatin_RNA = length(
    cis_rna_ids
  ),

  Final_Triple_Matched = length(
    matched_ids
  ),

  OralCavity_Final = nrow(
    oral_final
  ),

  Ultra7_Detected = length(
    ultra_present
  ),

  Ultra7_Total = length(
    ultra7
  ),

  Primary4_Detected = length(
    primary_present
  ),

  Primary4_Total = length(
    primary4
  ),

  Priority52_Detected = length(
    priority_present
  ),

  Priority52_ValidForScore = length(
    priority_valid
  ),

  Priority52_Total = length(
    priority52
  ),

  stringsAsFactors = FALSE
)

write.table(
  final,
  file = file.path(
    out_dir,
    "GDSC1_HNSCC_Cisplatin_Integrated_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  oral_final,
  file = file.path(
    out_dir,
    "GDSC1_OralCavity_Cisplatin_Integrated_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  qc_summary,
  file = file.path(
    out_dir,
    "GDSC1_QC_Summary_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  primary_result,
  file = file.path(
    out_dir,
    "GDSC1_PRIMARY_Ultra7_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  secondary_primary4,
  file = file.path(
    out_dir,
    "GDSC1_SECONDARY_Primary4_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ultra7_single,
  file = file.path(
    out_dir,
    "GDSC1_SENSITIVITY_Ultra7_SingleGenes_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  score_sensitivity,
  file = file.path(
    out_dir,
    "GDSC1_SENSITIVITY_Scores_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  oral_results,
  file = file.path(
    out_dir,
    "GDSC1_SENSITIVITY_OralCavity_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

step7b_results <- list(
  QC = qc_summary,
  Primary = primary_result,
  Secondary_Primary4 = secondary_primary4,
  Ultra7_SingleGenes = ultra7_single,
  Score_Sensitivity = score_sensitivity,
  Oral_Sensitivity = oral_results,
  HNSCC_Integrated = final,
  Oral_Integrated = oral_final,
  Frozen_Ultra7 = ultra7,
  Frozen_Primary4 = primary4,
  Frozen_Priority52 = priority52,
  Priority52_Used = priority_valid,
  Cisplatin_DRUG_ID = 1005
)

saveRDS(
  step7b_results,
  file = file.path(
    out_dir,
    "Part3_Step7B_GDSC1_FINAL_V3.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step7B_sessionInfo_FINAL_V3.txt"
  )
)

cat("\n")
cat("====================================================\n")
cat("PART3 STEP 7B GDSC1 FINAL V3 COMPLETE\n")
cat("====================================================\n\n")

cat("QC:\n")
print(qc_summary)

cat("\nPRIMARY ENDPOINT:\n")
print(primary_result)

cat("\nSECONDARY PRIMARY4:\n")
print(secondary_primary4)

cat("\nULTRA7 SINGLE-GENE SENSITIVITY:\n")
print(ultra7_single)

cat("\nSCORE SENSITIVITY:\n")
print(score_sensitivity)

cat("\nORAL CAVITY SENSITIVITY:\n")
print(oral_results)
