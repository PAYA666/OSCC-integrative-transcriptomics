project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

if (
  !exists("pairwise_results") ||
  !exists("kw_results") ||
  !exists("packyear_results") ||
  !exists("status_summary")
) {

  sens <- readRDS(
    file.path(
      proj_dir,
      "Part2_Smoking_ExposureDose_Sensitivity_FINAL.rds"
    )
  )

  status_summary   <- sens$status_summary
  kw_results       <- sens$kw_results
  pairwise_results <- sens$pairwise_results
  packyear_results <- sens$packyear_results
}

key_genes <- c(
  "CYP1A1",
  "CES1",
  "SFRP4",
  "SFRP2",
  "PDGFRA",
  "OGN",
  "FNDC1",
  "TK1",
  "MYBL2"
)

get_pair <- function(comp, prefix) {

  x <- pairwise_results[
    pairwise_results$Comparison == comp,
    c(
      "gene_name",
      "Median_Difference",
      "P",
      "FDR"
    )
  ]

  colnames(x)[-1] <- paste0(
    prefix,
    c(
      "_MedianDiff",
      "_P",
      "_FDR"
    )
  )

  x
}

CN <- get_pair(
  "Current_vs_Never",
  "Current_vs_Never"
)

FN <- get_pair(
  "Former_vs_Never",
  "Former_vs_Never"
)

CF <- get_pair(
  "Current_vs_Former",
  "Current_vs_Former"
)

summary_table <- merge(
  CN,
  FN,
  by = "gene_name",
  all = TRUE
)

summary_table <- merge(
  summary_table,
  CF,
  by = "gene_name",
  all = TRUE
)

kw2 <- kw_results[
  ,
  c(
    "gene_name",
    "P",
    "FDR"
  )
]

colnames(kw2) <- c(
  "gene_name",
  "Smoking3_Global_P",
  "Smoking3_Global_FDR"
)

summary_table <- merge(
  summary_table,
  kw2,
  by = "gene_name",
  all = TRUE
)

pack2 <- packyear_results[
  ,
  c(
    "gene_name",
    "N",
    "Spearman_Rho",
    "P",
    "FDR"
  )
]

colnames(pack2) <- c(
  "gene_name",
  "PackYears_N",
  "PackYears_Rho",
  "PackYears_P",
  "PackYears_FDR"
)

summary_table <- merge(
  summary_table,
  pack2,
  by = "gene_name",
  all = TRUE
)

median_wide <- reshape(
  status_summary[
    ,
    c(
      "gene_name",
      "Smoking_status",
      "Median"
    )
  ],
  idvar = "gene_name",
  timevar = "Smoking_status",
  direction = "wide"
)

colnames(median_wide) <- sub(
  "^Median\\.",
  "",
  colnames(median_wide)
)

colnames(median_wide)[
  colnames(median_wide) == "Never"
] <- "Never_Median"

colnames(median_wide)[
  colnames(median_wide) == "Former"
] <- "Former_Median"

colnames(median_wide)[
  colnames(median_wide) == "Current"
] <- "Current_Median"

summary_table <- merge(
  summary_table,
  median_wide,
  by = "gene_name",
  all = TRUE
)

summary_table$Current_vs_Never_sig <-
  summary_table$Current_vs_Never_FDR < 0.05

summary_table$Former_vs_Never_sig <-
  summary_table$Former_vs_Never_FDR < 0.05

summary_table$Current_vs_Former_sig <-
  summary_table$Current_vs_Former_FDR < 0.05

summary_table$PackYears_sig <-
  summary_table$PackYears_FDR < 0.05

summary_table$Pattern <- "No_clear_3group_pattern"

for (i in seq_len(nrow(summary_table))) {

  cn <- isTRUE(
    summary_table$Current_vs_Never_sig[i]
  )

  fn <- isTRUE(
    summary_table$Former_vs_Never_sig[i]
  )

  cf <- isTRUE(
    summary_table$Current_vs_Former_sig[i]
  )

  if (cn && !fn && cf) {

    summary_table$Pattern[i] <-
      "Current_specific"

  } else if (cn && !fn && !cf) {

    summary_table$Pattern[i] <-
      "Current_vs_Never_only"

  } else if (cn && fn && !cf) {

    summary_table$Pattern[i] <-
      "Persistent_Ever_association"

  } else if (cn && fn && cf) {

    summary_table$Pattern[i] <-
      "Current_stronger_than_Former"

  } else if (!cn && fn) {

    summary_table$Pattern[i] <-
      "Former_association_only"
  }
}

summary_table$gene_name <- factor(
  summary_table$gene_name,
  levels = key_genes
)

summary_table <- summary_table[
  order(summary_table$gene_name),
]

summary_table$gene_name <- as.character(
  summary_table$gene_name
)

write.table(
  summary_table,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_StatusPattern_FINAL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  summary_table,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_StatusPattern_FINAL.rds"
  )
)

cat("\n=====================================\n")
cat("Step 14B finished\n")
cat("=====================================\n")

print(
  summary_table[
    ,
    c(
      "gene_name",
      "Never_Median",
      "Former_Median",
      "Current_Median",
      "Current_vs_Never_FDR",
      "Former_vs_Never_FDR",
      "Current_vs_Former_FDR",
      "PackYears_Rho",
      "PackYears_FDR",
      "Pattern"
    )
  ]
)
