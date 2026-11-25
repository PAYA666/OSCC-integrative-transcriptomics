rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

root <- project_root

part3_dir <- file.path(
  root,
  "Part3_Cisplatin"
)

step3_file <- file.path(
  part3_dir,
  "Step3_Gene_Convergence",
  "Part3_Step3_GeneConvergence_COMPLETE.rds"
)

step5_file <- file.path(
  part3_dir,
  "Step5_TCGA_WGCNA_Robustness",
  "Part3_Step5_TCGA_WGCNA_Robustness_COMPLETE.rds"
)

required_files <- c(
  step3_file,
  step5_file
)

if (!all(file.exists(required_files))) {

  cat("\nMissing files:\n")

  print(
    required_files[
      !file.exists(required_files)
    ]
  )

  stop(
    "Required Step 3 / Step 5 file missing."
  )
}

out_dir <- file.path(
  part3_dir,
  "Step6_Candidate_Consolidation"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

step3 <- readRDS(
  step3_file
)

step5 <- readRDS(
  step5_file
)

required_step3 <- c(
  "priority_pool",
  "ultra_strict",
  "strict_both_all"
)

missing_step3 <- setdiff(
  required_step3,
  names(step3)
)

if (length(missing_step3) > 0) {

  stop(
    paste0(
      "Missing Step 3 object(s): ",
      paste(
        missing_step3,
        collapse = ", "
      )
    )
  )
}

priority_table <- as.data.frame(
  step3$priority_pool
)

ultra_table <- as.data.frame(
  step3$ultra_strict
)

strict_table <- as.data.frame(
  step3$strict_both_all
)

required_step5 <- c(
  "candidate_network",
  "priority_depletion",
  "strict_depletion"
)

missing_step5 <- setdiff(
  required_step5,
  names(step5)
)

if (length(missing_step5) > 0) {

  stop(
    paste0(
      "Missing Step 5 object(s): ",
      paste(
        missing_step5,
        collapse = ", "
      )
    )
  )
}

candidate_network <- as.data.frame(
  step5$candidate_network
)

priority_depletion <- as.data.frame(
  step5$priority_depletion
)

strict_depletion <- as.data.frame(
  step5$strict_depletion
)

for (obj_name in c(
  "priority_table",
  "ultra_table",
  "strict_table",
  "candidate_network"
)) {

  obj <- get(obj_name)

  if (!("Gene" %in% colnames(obj))) {

    stop(
      paste0(
        obj_name,
        " does not contain a Gene column."
      )
    )
  }
}

cat("\n========== INPUT COUNTS ==========\n")

cat(
  "Priority =",
  length(unique(priority_table$Gene)),
  "\n"
)

cat(
  "Ultra =",
  length(unique(ultra_table$Gene)),
  "\n"
)

cat(
  "StrictBoth =",
  length(unique(strict_table$Gene)),
  "\n"
)

evidence <- candidate_network

evidence <- evidence[
  !duplicated(evidence$Gene),
  ,
  drop = FALSE
]

fields_to_recover <- c(
  "logFC_197",
  "logFC_111",
  "Meta_Z",
  "Meta_FDR",
  "Pathway_Tags"
)

for (nm in fields_to_recover) {

  if (
    !(nm %in% colnames(evidence)) &&
    nm %in% colnames(priority_table)
  ) {

    evidence[, nm] <- priority_table[
      match(
        evidence$Gene,
        priority_table$Gene
      ),
      nm
    ]
  }
}

evidence$In_Priority52 <- evidence$Gene %in%
  priority_table$Gene

evidence$In_StrictBoth97 <- evidence$Gene %in%
  strict_table$Gene

evidence$In_Ultra7 <- evidence$Gene %in%
  ultra_table$Gene

if (!("module_color" %in% colnames(evidence))) {
  stop(
    "candidate_network does not contain module_color."
  )
}

evidence$Network_Interpretable <-
  !is.na(evidence$module_color) &
  evidence$module_color != "" &
  evidence$module_color != "grey"

if (!("Pathway_Tags" %in% colnames(evidence))) {

  evidence$Pathway_Tags <- NA_character_

}

count_pathways <- function(x) {

  if (
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(0L)
  }

  tags <- strsplit(
    x,
    ";",
    fixed = TRUE
  )[[1]]

  tags <- trimws(tags)

  tags <- tags[
    tags != ""
  ]

  length(
    unique(tags)
  )
}

evidence$Core_Pathway_Breadth <- vapply(
  evidence$Pathway_Tags,
  count_pathways,
  integer(1)
)

required_dep_cols <- c(
  "Module",
  "Depleted_Rho",
  "Depleted_FDR"
)

if (!all(
  required_dep_cols %in%
    colnames(priority_depletion)
)) {

  stop(
    "priority_depletion is missing required columns."
  )
}

if (!all(
  required_dep_cols %in%
    colnames(strict_depletion)
)) {

  stop(
    "strict_depletion is missing required columns."
  )
}

evidence$Priority_ModuleDepleted_Rho <-
  priority_depletion$Depleted_Rho[
    match(
      evidence$module_color,
      priority_depletion$Module
    )
  ]

evidence$Priority_ModuleDepleted_FDR <-
  priority_depletion$Depleted_FDR[
    match(
      evidence$module_color,
      priority_depletion$Module
    )
  ]

evidence$Strict_ModuleDepleted_Rho <-
  strict_depletion$Depleted_Rho[
    match(
      evidence$module_color,
      strict_depletion$Module
    )
  ]

evidence$Strict_ModuleDepleted_FDR <-
  strict_depletion$Depleted_FDR[
    match(
      evidence$module_color,
      strict_depletion$Module
    )
  ]

evidence$Priority_Module_Robust <-
  evidence$Network_Interpretable &
  is.finite(
    evidence$Priority_ModuleDepleted_FDR
  ) &
  evidence$Priority_ModuleDepleted_FDR < 0.05

evidence$Strict_Module_Robust <-
  evidence$Network_Interpretable &
  is.finite(
    evidence$Strict_ModuleDepleted_FDR
  ) &
  evidence$Strict_ModuleDepleted_FDR < 0.05

evidence$Candidate_Tier <- "Tier3_Priority"

evidence$Candidate_Tier[
  evidence$In_Ultra7 &
  !evidence$Network_Interpretable
] <- "Tier2_PathwayAnchored_Ultra"

evidence$Candidate_Tier[
  evidence$In_Ultra7 &
  evidence$Network_Interpretable
] <- "Tier1_NetworkAnchored_Ultra"

tier1 <- evidence[
  evidence$Candidate_Tier ==
    "Tier1_NetworkAnchored_Ultra",
  ,
  drop = FALSE
]

tier2 <- evidence[
  evidence$Candidate_Tier ==
    "Tier2_PathwayAnchored_Ultra",
  ,
  drop = FALSE
]

if (!("Meta_FDR" %in% colnames(tier1))) {
  stop(
    "Meta_FDR is required for frozen lead-candidate rule."
  )
}

if (!("LOO_Rho" %in% colnames(tier1))) {
  stop(
    "LOO_Rho is required for frozen lead-candidate rule."
  )
}

tier1_order <- order(
  -tier1$Core_Pathway_Breadth,
  tier1$Meta_FDR,
  -abs(tier1$LOO_Rho),
  na.last = TRUE
)

tier1 <- tier1[
  tier1_order,
  ,
  drop = FALSE
]

if (nrow(tier1) > 0) {

  lead_gene <- tier1$Gene[1]

} else {

  lead_gene <- NA_character_

}

evidence$Lead_Candidate <-
  evidence$Gene == lead_gene

tier_rank <- c(
  "Tier1_NetworkAnchored_Ultra" = 1,
  "Tier2_PathwayAnchored_Ultra" = 2,
  "Tier3_Priority" = 3
)

evidence$Tier_Rank <- unname(
  tier_rank[
    evidence$Candidate_Tier
  ]
)

evidence <- evidence[
  order(
    evidence$Tier_Rank,
    -evidence$Core_Pathway_Breadth,
    evidence$Meta_FDR,
    -abs(evidence$LOO_Rho),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

ultra_evidence <- evidence[
  evidence$In_Ultra7,
  ,
  drop = FALSE
]

publication_cols <- c(
  "Gene",
  "Candidate_Tier",
  "Lead_Candidate",
  "logFC_197",
  "logFC_111",
  "Meta_Z",
  "Meta_FDR",
  "Pathway_Tags",
  "Core_Pathway_Breadth",
  "module_color",
  "Assigned_kME",
  "LOO_Rho",
  "LOO_FDR",
  "Priority_ModuleDepleted_Rho",
  "Priority_ModuleDepleted_FDR",
  "Strict_ModuleDepleted_Rho",
  "Strict_ModuleDepleted_FDR"
)

publication_cols <- publication_cols[
  publication_cols %in%
    colnames(ultra_evidence)
]

ultra_publication <- ultra_evidence[
  ,
  publication_cols,
  drop = FALSE
]

primary_candidates <- tier1$Gene

secondary_candidates <- tier2$Gene

candidate_manifest <- data.frame(
  Category = c(
    "Lead_candidate",
    "Primary_candidate_panel",
    "Secondary_ultra_panel",
    "Discovery_basis",
    "TCGA_role",
    "Excluded_from_selection",
    "Future_external_validation_rule"
  ),

  Frozen_Definition = c(
    lead_gene,

    paste(
      primary_candidates,
      collapse = ";"
    ),

    paste(
      secondary_candidates,
      collapse = ";"
    ),

    paste0(
      "GSE197561 + GSE111585 independent acquired-resistance analyses; ",
      "Hallmark pathway convergence; signed gene-level convergence"
    ),

    paste0(
      "Biological/network projection only; ",
      "not clinical cisplatin-response validation"
    ),

    paste0(
      "GSE117872 and E-MTAB-15148 are supplementary/exploratory ",
      "and do not affect candidate ranking"
    ),

    paste0(
      "Any future external dataset must test the frozen Ultra7 / ",
      "primary candidate panel without re-ranking genes post hoc"
    )
  ),

  stringsAsFactors = FALSE
)

freeze_rules <- c(
  "PART 3 CANDIDATE FREEZE",
  "",
  "Candidate discovery and ranking use Step 1-5 evidence only.",
  "",
  "Tier 1 = Ultra-strict gene located in a non-grey frozen WGCNA module.",
  "Tier 2 = Ultra-strict gene located in grey/unassigned WGCNA.",
  "Tier 3 = remaining pathway-informed Priority52 genes.",
  "",
  "No arbitrary kME, LOO correlation, or weighted composite cutoff is introduced.",
  "",
  "Lead candidate rule within Tier 1:",
  "1. Maximum frozen core-pathway breadth.",
  "2. If tied, smaller Meta_FDR.",
  "3. If still tied, larger absolute leave-one-gene-out correlation.",
  "",
  paste0(
    "Frozen lead candidate: ",
    lead_gene
  ),
  "",
  paste0(
    "Frozen primary candidate panel: ",
    paste(
      primary_candidates,
      collapse = ", "
    )
  ),
  "",
  paste0(
    "Frozen secondary ultra-strict panel: ",
    paste(
      secondary_candidates,
      collapse = ", "
    )
  ),
  "",
  "GSE117872 and E-MTAB-15148 do not participate in candidate selection.",
  "",
  "TCGA/WGCNA findings are projection/network-context evidence, not clinical chemoresistance evidence."
)

writeLines(
  freeze_rules,
  con = file.path(
    out_dir,
    "Part3_Candidate_Freeze_Rules.txt"
  )
)

write.table(
  evidence,
  file = file.path(
    out_dir,
    "Part3_Priority52_Frozen_Evidence_Matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ultra_publication,
  file = file.path(
    out_dir,
    "Part3_Ultra7_Frozen_Evidence_Matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  candidate_manifest,
  file = file.path(
    out_dir,
    "Part3_Frozen_Candidate_Manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

step6_consolidation <- list(
  full_priority_evidence = evidence,
  ultra7_evidence = ultra_evidence,
  ultra7_publication = ultra_publication,
  tier1_primary = tier1,
  tier2_secondary = tier2,
  lead_gene = lead_gene,
  primary_candidates = primary_candidates,
  secondary_candidates = secondary_candidates,
  candidate_manifest = candidate_manifest,
  freeze_rules = freeze_rules
)

saveRDS(
  step6_consolidation,
  file = file.path(
    out_dir,
    "Part3_Step6_CandidateConsolidation_FROZEN.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step6_sessionInfo.txt"
  )
)

cat("\n")
cat("####################################################\n")
cat("STEP 6 CANDIDATE CONSOLIDATION FINISHED\n")
cat("####################################################\n\n")

cat(
  "Lead candidate:",
  lead_gene,
  "\n"
)

cat(
  "Primary candidate panel:",
  paste(
    primary_candidates,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Secondary ultra panel:",
  paste(
    secondary_candidates,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Output directory:",
  out_dir,
  "\n"
)
