rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

part3_dir <- file.path(project_root, "Part3_Cisplatin")

step6_file <- file.path(
  part3_dir,
  "Step6_Candidate_Consolidation",
  "Part3_Step6_CandidateConsolidation_FROZEN.rds"
)

step7c_file <- file.path(
  part3_dir,
  "GDSC",
  "Step7C_GDSC1_Result_Freeze",
  "Part3_Step7C_GDSC1_ResultFreeze_COMPLETE.rds"
)

out_dir <- file.path(
  part3_dir,
  "Step8_Final_Evidence_Integration"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(step6_file)) {
  stop("Step6 frozen RDS not found.")
}

if (!file.exists(step7c_file)) {
  stop("Step7C GDSC1 frozen RDS not found.")
}

step6 <- readRDS(step6_file)

step7c <- readRDS(step7c_file)

lead_gene <- "CCN1"

primary4 <- c(
  "CCN1",
  "CCL2",
  "ZBTB10",
  "DAB2"
)

secondary3 <- c(
  "JUN",
  "KYNU",
  "CCL20"
)

ultra7 <- c(
  "CCN1",
  "JUN",
  "CCL2",
  "ZBTB10",
  "CCL20",
  "DAB2",
  "KYNU"
)

final_candidate_order <- c(
  "CCN1",
  "CCL2",
  "ZBTB10",
  "DAB2",
  "JUN",
  "KYNU",
  "CCL20"
)

if (!"full_priority_evidence" %in% names(step6)) {
  stop("Step6 object lacks full_priority_evidence.")
}

evidence <- as.data.frame(
  step6[["full_priority_evidence"]]
)

if (!"Gene" %in% colnames(evidence)) {
  stop("Step6 evidence lacks Gene column.")
}

evidence$Gene <- as.character(
  evidence$Gene
)

ultra_evidence <- evidence[
  evidence$Gene %in% ultra7,
  ,
  drop = FALSE
]

ultra_evidence <- ultra_evidence[
  match(
    final_candidate_order,
    ultra_evidence$Gene
  ),
  ,
  drop = FALSE
]

if (
  nrow(ultra_evidence) != 7 ||
  any(is.na(ultra_evidence$Gene))
) {
  stop("Could not recover all seven frozen Ultra7 genes from Step6.")
}

ultra_evidence$Frozen_Class <- ifelse(
  ultra_evidence$Gene == lead_gene,
  "Lead",
  ifelse(
    ultra_evidence$Gene %in% primary4,
    "Primary",
    "Secondary"
  )
)

ultra_evidence$Frozen_Order <- match(
  ultra_evidence$Gene,
  final_candidate_order
)

if (!"Ultra7_SingleGenes" %in% names(step7c)) {
  stop("Step7C object lacks Ultra7_SingleGenes.")
}

gdsc_gene <- as.data.frame(
  step7c[["Ultra7_SingleGenes"]]
)

required_gene_cols <- c(
  "Gene",
  "N",
  "Spearman_Rho",
  "P_value",
  "Expected_Direction",
  "Direction_Supported",
  "FDR_BH7"
)

if (!all(required_gene_cols %in% colnames(gdsc_gene))) {
  stop("Unexpected GDSC1 Ultra7 single-gene result structure.")
}

gdsc_gene <- gdsc_gene[
  gdsc_gene$Gene %in% ultra7,
  required_gene_cols,
  drop = FALSE
]

colnames(gdsc_gene) <- c(
  "Gene",
  "GDSC1_N",
  "GDSC1_Rho",
  "GDSC1_P",
  "GDSC1_ExpectedDirection",
  "GDSC1_DirectionSupported",
  "GDSC1_FDR_BH7"
)

if (!"Primary4" %in% names(step7c)) {
  stop("Step7C object lacks Primary4 results.")
}

gdsc_primary4 <- as.data.frame(
  step7c[["Primary4"]]
)

required_primary4_cols <- c(
  "Gene",
  "FDR_BH4"
)

if (!all(required_primary4_cols %in% colnames(gdsc_primary4))) {
  stop("Unexpected GDSC1 Primary4 result structure.")
}

gdsc_primary4 <- gdsc_primary4[
  gdsc_primary4$Gene %in% primary4,
  c(
    "Gene",
    "FDR_BH4"
  ),
  drop = FALSE
]

colnames(gdsc_primary4) <- c(
  "Gene",
  "GDSC1_FDR_BH4"
)

final_evidence <- merge(
  ultra_evidence,
  gdsc_gene,
  by = "Gene",
  all.x = TRUE,
  sort = FALSE
)

final_evidence <- merge(
  final_evidence,
  gdsc_primary4,
  by = "Gene",
  all.x = TRUE,
  sort = FALSE
)

final_evidence <- final_evidence[
  order(final_evidence$Frozen_Order),
  ,
  drop = FALSE
]

rownames(final_evidence) <- NULL

final_evidence$GDSC1_Evidence <- NA_character_

for (i in seq_len(nrow(final_evidence))) {

  gene <- final_evidence$Gene[i]

  rho <- final_evidence$GDSC1_Rho[i]

  pval <- final_evidence$GDSC1_P[i]

  direction_ok <- final_evidence$GDSC1_DirectionSupported[i]

  if (gene %in% primary4) {

    fdr4 <- final_evidence$GDSC1_FDR_BH4[i]

    if (
      !is.na(fdr4) &&
      fdr4 < 0.05 &&
      isTRUE(direction_ok)
    ) {

      final_evidence$GDSC1_Evidence[i] <-
        "BH4-significant directional support"

    } else if (
      !is.na(pval) &&
      pval < 0.05 &&
      isTRUE(direction_ok)
    ) {

      final_evidence$GDSC1_Evidence[i] <-
        "Nominal directional support; not BH4-significant"

    } else if (isTRUE(direction_ok)) {

      final_evidence$GDSC1_Evidence[i] <-
        "Direction supported; not significant"

    } else {

      final_evidence$GDSC1_Evidence[i] <-
        "Expected direction not supported"
    }

  } else {

    fdr7 <- final_evidence$GDSC1_FDR_BH7[i]

    if (
      !is.na(fdr7) &&
      fdr7 < 0.05 &&
      isTRUE(direction_ok)
    ) {

      final_evidence$GDSC1_Evidence[i] <-
        "BH7-significant directional support"

    } else if (
      !is.na(pval) &&
      pval < 0.05 &&
      isTRUE(direction_ok)
    ) {

      final_evidence$GDSC1_Evidence[i] <-
        "Nominal directional support"

    } else if (isTRUE(direction_ok)) {

      final_evidence$GDSC1_Evidence[i] <-
        "Direction supported; not significant"

    } else {

      final_evidence$GDSC1_Evidence[i] <-
        "Expected direction not supported"
    }
  }
}

final_evidence$Candidate_Status_After_GDSC1 <- ifelse(
  final_evidence$Gene == lead_gene,
  "Lead remains frozen",
  ifelse(
    final_evidence$Gene %in% primary4,
    "Primary candidate remains frozen",
    "Secondary Ultra7 candidate remains frozen"
  )
)

if (!"Primary" %in% names(step7c)) {
  stop("Step7C object lacks Primary endpoint.")
}

gdsc_primary <- as.data.frame(
  step7c[["Primary"]]
)

required_primary_cols <- c(
  "N",
  "Spearman_Rho",
  "P_value"
)

if (!all(required_primary_cols %in% colnames(gdsc_primary))) {
  stop("Unexpected GDSC1 primary endpoint structure.")
}

primary_summary <- data.frame(
  Endpoint =
    "Frozen Ultra7 score vs GDSC1 Cisplatin LN_IC50",

  Cohort =
    "HNSCC",

  N =
    gdsc_primary$N[1],

  Spearman_Rho =
    gdsc_primary$Spearman_Rho[1],

  P_value =
    gdsc_primary$P_value[1],

  Prespecified_Direction =
    "Positive",

  Result = ifelse(
    gdsc_primary$P_value[1] < 0.05 &&
      gdsc_primary$Spearman_Rho[1] > 0,
    "Statistically supported",
    ifelse(
      gdsc_primary$Spearman_Rho[1] > 0,
      "Positive direction but not statistically significant",
      "Expected direction not supported"
    )
  ),

  stringsAsFactors = FALSE
)

evidence_hierarchy <- data.frame(
  Evidence_Level = c(
    "1",
    "2",
    "3",
    "4"
  ),

  Evidence_Source = c(
    "Two independent acquired cisplatin-resistance bulk datasets",
    "Cross-dataset gene and pathway convergence",
    "TCGA OSCC and frozen WGCNA biological projection",
    "GDSC1 HNSCC pharmacogenomic assessment"
  ),

  Interpretation = c(
    paste0(
      "Primary discovery and independent replication evidence ",
      "for experimentally acquired cisplatin resistance."
    ),

    paste0(
      "Defines the frozen resistance-associated candidate ",
      "architecture, including the Ultra7 set."
    ),

    paste0(
      "Shows tumor-context and network coordination; ",
      "does not represent clinical cisplatin resistance."
    ),

    paste0(
      "Tests whether the frozen acquired-resistance program ",
      "generalizes to baseline cross-cell-line intrinsic ",
      "cisplatin sensitivity."
    )
  ),

  Final_Status = c(
    "Core evidence",
    "Core evidence",
    "Supporting biological context",
    "Orthogonal assessment; Ultra7 primary endpoint not significant"
  ),

  stringsAsFactors = FALSE
)

final_conclusions <- data.frame(
  Item = c(
    "Core_resistance_program",
    "Candidate_freeze",
    "Lead_gene",
    "GDSC1_primary",
    "GDSC1_DAB2",
    "GDSC1_interpretation",
    "Final_boundary"
  ),

  Conclusion = c(
    paste0(
      "The principal cross-model acquired cisplatin-resistance ",
      "signal is a conserved stress-adaptation program dominated ",
      "by NF-kB/hypoxia-associated signaling, with EMT/TGF-beta-related ",
      "plasticity as supportive biology."
    ),

    paste0(
      "The frozen Ultra7 candidate set is CCN1, JUN, CCL2, ZBTB10, ",
      "CCL20, DAB2 and KYNU; the primary candidate panel is CCN1, ",
      "CCL2, ZBTB10 and DAB2."
    ),

    paste0(
      "CCN1 remains the prespecified lead candidate based on Step1-6 ",
      "evidence. No later external dataset is used to rerank the ",
      "frozen candidates."
    ),

    paste0(
      "In GDSC1 HNSCC models, the frozen Ultra7 score showed a weak ",
      "positive but nonsignificant association with Cisplatin LN_IC50 ",
      "(rho=",
      sprintf("%.3f", gdsc_primary$Spearman_Rho[1]),
      ", P=",
      signif(gdsc_primary$P_value[1], 3),
      ")."
    ),

    paste0(
      "DAB2 showed the strongest single-gene pharmacogenomic signal ",
      "among the frozen primary candidates, but only at the nominal level ",
      "and did not remain significant after BH correction across Primary4."
    ),

    paste0(
      "GDSC1 therefore does not statistically validate Ultra7 as a marker ",
      "of cross-cell-line intrinsic cisplatin sensitivity. This does not ",
      "invalidate the acquired-resistance signal because acquired adaptation ",
      "and baseline intrinsic sensitivity are biologically distinct phenotypes."
    ),

    paste0(
      "The computational candidate-selection main line is considered complete ",
      "after the prespecified GDSC1 assessment. No post-hoc candidate ",
      "optimization is performed."
    )
  ),

  stringsAsFactors = FALSE
)

label_correction <- data.frame(
  Incorrect = "Primary4_CCNI",
  Correct = "Primary4_CCN1",
  Note = "Label typo only; no numerical result was affected.",
  stringsAsFactors = FALSE
)

if (nrow(final_conclusions) != 7) {
  stop("Final conclusion table should contain exactly 7 rows.")
}

if (nrow(final_evidence) != 7) {
  stop("Final Ultra7 evidence table should contain exactly 7 genes.")
}

if (!identical(
  final_evidence$Gene,
  final_candidate_order
)) {
  stop("Final candidate order is incorrect.")
}

if (lead_gene != "CCN1") {
  stop("Frozen lead gene unexpectedly changed.")
}

write.table(
  final_evidence,
  file = file.path(
    out_dir,
    "Part3_FINAL_Ultra7_Evidence_Table_V2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  primary_summary,
  file = file.path(
    out_dir,
    "Part3_FINAL_GDSC1_PrimaryEndpoint_V2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  evidence_hierarchy,
  file = file.path(
    out_dir,
    "Part3_FINAL_Evidence_Hierarchy_V2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  final_conclusions,
  file = file.path(
    out_dir,
    "Part3_FINAL_Conclusions_V2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  label_correction,
  file = file.path(
    out_dir,
    "Part3_Label_Correction_V2.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_lines <- c(
  "PART 3: OSCC CISPLATIN RESISTANCE - FINAL FROZEN SUMMARY",
  "",
  "CORE ACQUIRED-RESISTANCE EVIDENCE",
  "- GSE197561: discovery dataset",
  "- GSE111585: independent replication dataset",
  "- Independent DEG and Hallmark analyses",
  "- Cross-dataset gene-level and pathway-level convergence",
  "",
  "FROZEN CANDIDATES",
  paste0(
    "- Lead: ",
    lead_gene
  ),
  paste0(
    "- Primary4: ",
    paste(
      primary4,
      collapse = ", "
    )
  ),
  paste0(
    "- Secondary Ultra7: ",
    paste(
      secondary3,
      collapse = ", "
    )
  ),
  paste0(
    "- Ultra7: ",
    paste(
      ultra7,
      collapse = ", "
    )
  ),
  "",
  "GDSC1 PHARMACOGENOMIC ASSESSMENT",
  paste0(
    "- HNSCC n = ",
    gdsc_primary$N[1]
  ),
  paste0(
    "- Ultra7 vs Cisplatin LN_IC50: rho = ",
    sprintf(
      "%.3f",
      gdsc_primary$Spearman_Rho[1]
    ),
    ", P = ",
    signif(
      gdsc_primary$P_value[1],
      3
    )
  ),
  "- Prespecified Ultra7 primary endpoint was not statistically significant.",
  "- DAB2 showed nominal positive association but not BH4-significant.",
  "",
  "INTERPRETIVE BOUNDARY",
  "- GDSC1 assesses baseline cross-cell-line intrinsic cisplatin sensitivity.",
  "- GDSC1 is not a direct replication of experimentally acquired resistance.",
  "- GDSC1 does not alter the frozen candidate hierarchy.",
  "- CCN1 remains the frozen lead candidate.",
  "- No post-hoc candidate optimization is performed."
)

writeLines(
  summary_lines,
  con = file.path(
    out_dir,
    "Part3_FINAL_Frozen_Summary_V2.txt"
  )
)

step8 <- list(
  Lead = lead_gene,
  Primary4 = primary4,
  Secondary3 = secondary3,
  Ultra7 = ultra7,
  Final_Evidence = final_evidence,
  GDSC1_Primary = primary_summary,
  Evidence_Hierarchy = evidence_hierarchy,
  Final_Conclusions = final_conclusions,
  Label_Correction = label_correction
)

saveRDS(
  step8,
  file = file.path(
    out_dir,
    "Part3_Step8_FINAL_EvidenceIntegration_FROZEN_V2.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step8_sessionInfo_V2.txt"
  )
)

cat("\n")
cat("====================================================\n")
cat("PART 3 STEP 8 FINAL V2 COMPLETE\n")
cat("====================================================\n\n")

cat(
  "Lead:",
  lead_gene,
  "\n"
)

cat(
  "Primary4:",
  paste(
    primary4,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Ultra7:",
  paste(
    ultra7,
    collapse = ", "
  ),
  "\n\n"
)

print(primary_summary)

cat("\nFinal Ultra7 evidence:\n")

print(
  final_evidence[
    ,
    c(
      "Gene",
      "Frozen_Class",
      "GDSC1_Rho",
      "GDSC1_P",
      "GDSC1_FDR_BH7",
      "GDSC1_FDR_BH4",
      "GDSC1_Evidence",
      "Candidate_Status_After_GDSC1"
    ),
    drop = FALSE
  ]
)

cat("\nFinal conclusions:\n")

print(final_conclusions)
