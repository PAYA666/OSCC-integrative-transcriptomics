rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)
options(timeout = 1200)

set.seed(20260815)

required_pkgs <- c(
  "fgsea",
  "ggplot2",
  "pheatmap"
)

pkg_check <- sapply(
  required_pkgs,
  requireNamespace,
  quietly = TRUE
)

print(pkg_check)

if (!all(pkg_check)) {
  stop(
    paste0(
      "Missing package(s): ",
      paste(
        names(pkg_check)[!pkg_check],
        collapse = ", "
      )
    )
  )
}

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(pheatmap)
})

project_dir <- file.path(project_root, "Part3_Cisplatin")

step1_file <- file.path(
  project_dir,
  "Step1_FINAL_V4_OFFLINE",
  "Part3_Step1_FINAL_V4_AllResults.rds"
)

if (!file.exists(step1_file)) {
  stop(
    paste0(
      "Step 1 RDS not found:\n",
      step1_file
    )
  )
}

gmt_name <- "h.all.v2026.1.Hs.symbols.gmt"

candidate_gmt <- c(
  file.path(
    project_dir,
    gmt_name
  ),
  file.path(
    project_root,
    gmt_name
  )
)

candidate_gmt <- candidate_gmt[
  file.exists(candidate_gmt)
]

if (length(candidate_gmt) == 0) {

  cat(
    "\nSearching the project root recursively for GMT...\n"
  )

  candidate_gmt <- list.files(
    path = project_root,
    pattern = "^h\\.all\\.v2026\\.1\\.Hs\\.symbols\\.gmt$",
    recursive = TRUE,
    full.names = TRUE
  )
}

if (length(candidate_gmt) == 0) {
  stop(
    paste0(
      "Cannot find ",
      gmt_name,
      " under the project root"
    )
  )
}

gmt_file <- candidate_gmt[1]

cat(
  "\nHallmark GMT:\n",
  gmt_file,
  "\n"
)

out_dir <- file.path(
  project_dir,
  "Step2_Hallmark_GSEA_V2"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

results <- readRDS(
  step1_file
)

stopifnot(
  all(
    c(
      "GSE197561",
      "GSE111585"
    ) %in%
      names(results)
  )
)

cat("\nStep 1 datasets:\n")
print(names(results))

hallmark <- fgsea::gmtPathways(
  gmt_file
)

cat(
  "\nHallmark gene sets:",
  length(hallmark),
  "\n"
)

print(
  head(
    names(hallmark),
    10
  )
)

make_rank <- function(
    deg_table,
    dataset
) {

  needed <- c(
    "Gene",
    "t"
  )

  if (!all(needed %in% colnames(deg_table))) {
    stop(
      paste0(
        dataset,
        ": Gene or t column missing."
      )
    )
  }

  x <- deg_table[
    ,
    c(
      "Gene",
      "t"
    )
  ]

  x <- x[
    !is.na(x$Gene) &
      x$Gene != "" &
      is.finite(x$t),
    ,
    drop = FALSE
  ]

  x <- x[
    order(
      -abs(x$t)
    ),
    ,
    drop = FALSE
  ]

  x <- x[
    !duplicated(x$Gene),
    ,
    drop = FALSE
  ]

  ranks <- x$t

  names(ranks) <- x$Gene

  ranks <- sort(
    ranks,
    decreasing = TRUE
  )

  n_duplicate_stats <- sum(
    duplicated(ranks)
  )

  pct_duplicate_stats <- 100 *
    n_duplicate_stats /
    length(ranks)

  cat(
    dataset,
    "| ranked genes =",
    length(ranks),
    "| repeated statistic values =",
    n_duplicate_stats,
    sprintf(
      "(%.3f%%)\n",
      pct_duplicate_stats
    )
  )

  return(
    list(
      ranks = ranks,
      n_duplicate_stats =
        n_duplicate_stats,
      pct_duplicate_stats =
        pct_duplicate_stats
    )
  )
}

fgsea_export_table <- function(x) {

  y <- as.data.frame(
    x,
    stringsAsFactors = FALSE
  )

  list_columns <- vapply(
    y,
    is.list,
    logical(1)
  )

  if (any(list_columns)) {

    for (nm in names(y)[list_columns]) {

      y[[nm]] <- vapply(
        y[[nm]],
        function(z) {

          if (length(z) == 0) {
            return("")
          }

          paste(
            as.character(z),
            collapse = ";"
          )
        },
        character(1)
      )
    }
  }

  return(y)
}

run_hallmark <- function(
    dataset,
    deg_table
) {

  cat("\n")
  cat("========================================\n")
  cat("Running Hallmark GSEA:", dataset, "\n")
  cat("========================================\n")

  rank_result <- make_rank(
    deg_table = deg_table,
    dataset = dataset
  )

  ranks <- rank_result$ranks

  warning_log <- character(0)

  fg <- withCallingHandlers(

    fgsea::fgseaMultilevel(
      pathways = hallmark,
      stats = ranks,
      minSize = 15,
      maxSize = 500,
      eps = 0
    ),

    warning = function(w) {

      warning_log <<- c(
        warning_log,
        conditionMessage(w)
      )

      invokeRestart(
        "muffleWarning"
      )
    }
  )

  fg <- as.data.frame(
    fg,
    stringsAsFactors = FALSE
  )

  fg <- fg[
    order(
      fg$padj,
      -abs(fg$NES)
    ),
    ,
    drop = FALSE
  ]

  fg$Direction <- ifelse(
    fg$NES > 0,
    "Resistant",
    "Parental"
  )

  fg$Significant <- ifelse(
    !is.na(fg$padj) &
      fg$padj < 0.05,
    "FDR<0.05",
    "NS"
  )

  fg_export <- fgsea_export_table(
    fg
  )

  write.table(
    fg_export,
    file = file.path(
      out_dir,
      paste0(
        dataset,
        "_Hallmark_GSEA_ALL.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  fg_sig <- fg[
    !is.na(fg$padj) &
      fg$padj < 0.05,
    ,
    drop = FALSE
  ]

  fg_sig_export <- fgsea_export_table(
    fg_sig
  )

  write.table(
    fg_sig_export,
    file = file.path(
      out_dir,
      paste0(
        dataset,
        "_Hallmark_GSEA_FDR005.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  warning_log <- unique(
    warning_log
  )

  if (length(warning_log) == 0) {

    writeLines(
      "No fgsea warnings.",
      con = file.path(
        out_dir,
        paste0(
          dataset,
          "_fgsea_warnings.txt"
        )
      )
    )

  } else {

    writeLines(
      warning_log,
      con = file.path(
        out_dir,
        paste0(
          dataset,
          "_fgsea_warnings.txt"
        )
      )
    )
  }

  valid_plot <- fg[
    is.finite(fg$NES) &
      !is.na(fg$padj),
    ,
    drop = FALSE
  ]

  valid_plot <- valid_plot[
    order(
      valid_plot$padj,
      -abs(valid_plot$NES)
    ),
    ,
    drop = FALSE
  ]

  plot_df <- head(
    valid_plot,
    20
  )

  plot_df$pathway <- factor(
    plot_df$pathway,
    levels = rev(
      plot_df$pathway
    )
  )

  p <- ggplot(
    plot_df,
    aes(
      x = NES,
      y = pathway
    )
  ) +
    geom_col() +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    theme_classic(
      base_size = 11
    ) +
    labs(
      title = paste0(
        dataset,
        " Hallmark GSEA"
      ),
      x = "Normalized Enrichment Score",
      y = NULL
    )

  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        dataset,
        "_Hallmark_Top20.pdf"
      )
    ),
    plot = p,
    width = 8,
    height = 7
  )

  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        dataset,
        "_Hallmark_Top20.png"
      )
    ),
    plot = p,
    width = 8,
    height = 7,
    dpi = 300
  )

  saveRDS(
    fg,
    file = file.path(
      out_dir,
      paste0(
        dataset,
        "_Hallmark_GSEA_FULL_OBJECT.rds"
      )
    )
  )

  return(
    list(
      ranks = ranks,
      rank_qc = rank_result,
      all = fg,
      sig = fg_sig,
      warnings = warning_log
    )
  )
}

gsea197 <- run_hallmark(
  dataset = "GSE197561",
  deg_table =
    results$GSE197561$deg_all
)

gsea111 <- run_hallmark(
  dataset = "GSE111585",
  deg_table =
    results$GSE111585$deg_all
)

x197 <- gsea197$all[
  ,
  c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size"
  )
]

colnames(x197) <- c(
  "Pathway",
  "NES_GSE197561",
  "P_GSE197561",
  "FDR_GSE197561",
  "Size_GSE197561"
)

x111 <- gsea111$all[
  ,
  c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size"
  )
]

colnames(x111) <- c(
  "Pathway",
  "NES_GSE111585",
  "P_GSE111585",
  "FDR_GSE111585",
  "Size_GSE111585"
)

concordance <- merge(
  x197,
  x111,
  by = "Pathway",
  all = TRUE
)

valid_nes <- is.finite(
  concordance$NES_GSE197561
) &
  is.finite(
    concordance$NES_GSE111585
  )

concordance$Same_direction <- FALSE

concordance$Same_direction[
  valid_nes
] <- sign(
  concordance$NES_GSE197561[
    valid_nes
  ]
) ==
  sign(
    concordance$NES_GSE111585[
      valid_nes
    ]
  )

sig197 <- !is.na(
  concordance$FDR_GSE197561
) &
  concordance$FDR_GSE197561 < 0.05

sig111 <- !is.na(
  concordance$FDR_GSE111585
) &
  concordance$FDR_GSE111585 < 0.05

concordance$N_FDR005 <-
  as.integer(sig197) +
  as.integer(sig111)

concordance$Direction <- "Unavailable"

concordance$Direction[
  valid_nes &
    concordance$Same_direction &
    concordance$NES_GSE197561 > 0
] <- "Resistant"

concordance$Direction[
  valid_nes &
    concordance$Same_direction &
    concordance$NES_GSE197561 < 0
] <- "Parental"

concordance$Direction[
  valid_nes &
    !concordance$Same_direction
] <- "Opposite"

concordance$Evidence_Class <-
  "Unavailable"

concordance$Evidence_Class[
  valid_nes &
    concordance$Same_direction &
    concordance$N_FDR005 == 2
] <- "Concordant_FDR_both"

concordance$Evidence_Class[
  valid_nes &
    concordance$Same_direction &
    concordance$N_FDR005 == 1
] <- "Concordant_FDR_one"

concordance$Evidence_Class[
  valid_nes &
    concordance$Same_direction &
    concordance$N_FDR005 == 0
] <- "Same_direction_NS"

concordance$Evidence_Class[
  valid_nes &
    !concordance$Same_direction
] <- "Opposite_direction"

concordance$Mean_NES <- rowMeans(
  concordance[
    ,
    c(
      "NES_GSE197561",
      "NES_GSE111585"
    )
  ],
  na.rm = TRUE
)

concordance$Min_FDR <- apply(
  concordance[
    ,
    c(
      "FDR_GSE197561",
      "FDR_GSE111585"
    )
  ],
  1,
  function(z) {

    z <- z[
      is.finite(z)
    ]

    if (length(z) == 0) {
      return(NA_real_)
    }

    min(z)
  }
)

evidence_order <- c(
  "Concordant_FDR_both",
  "Concordant_FDR_one",
  "Same_direction_NS",
  "Opposite_direction",
  "Unavailable"
)

concordance$Evidence_Class <- factor(
  concordance$Evidence_Class,
  levels = evidence_order
)

concordance <- concordance[
  order(
    concordance$Evidence_Class,
    concordance$Min_FDR,
    -abs(concordance$Mean_NES),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

concordance$Evidence_Class <-
  as.character(
    concordance$Evidence_Class
  )

write.table(
  concordance,
  file = file.path(
    out_dir,
    "Part3_Hallmark_Pathway_Concordance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

strict_concordant <- concordance[
  concordance$Evidence_Class ==
    "Concordant_FDR_both",
  ,
  drop = FALSE
]

write.table(
  strict_concordant,
  file = file.path(
    out_dir,
    "Part3_Hallmark_Strict_Concordant_FDR_Both.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

supported_concordant <- concordance[
  concordance$Evidence_Class %in%
    c(
      "Concordant_FDR_both",
      "Concordant_FDR_one"
    ),
  ,
  drop = FALSE
]

write.table(
  supported_concordant,
  file = file.path(
    out_dir,
    "Part3_Hallmark_Supported_Concordant.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

key_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MTORC1_SIGNALING"
)

key_table <- concordance[
  concordance$Pathway %in%
    key_pathways,
  ,
  drop = FALSE
]

key_order <- match(
  key_table$Pathway,
  key_pathways
)

key_table <- key_table[
  order(key_order),
  ,
  drop = FALSE
]

write.table(
  key_table,
  file = file.path(
    out_dir,
    "Part3_Hallmark_PreSpecified_KeyPathways.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_concordance <- concordance[
  is.finite(
    concordance$NES_GSE197561
  ) &
    is.finite(
      concordance$NES_GSE111585
    ),
  ,
  drop = FALSE
]

label_df <- plot_concordance[
  plot_concordance$Evidence_Class %in%
    c(
      "Concordant_FDR_both",
      "Concordant_FDR_one"
    ),
  ,
  drop = FALSE
]

if (nrow(label_df) > 15) {

  label_df <- label_df[
    order(
      label_df$Min_FDR,
      -abs(label_df$Mean_NES)
    ),
    ,
    drop = FALSE
  ]

  label_df <- head(
    label_df,
    15
  )
}

p_scatter <- ggplot(
  plot_concordance,
  aes(
    x = NES_GSE197561,
    y = NES_GSE111585
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_point(
    size = 2.8
  ) +
  geom_text(
    data = label_df,
    aes(
      label = Pathway
    ),
    size = 2.5,
    check_overlap = TRUE,
    vjust = -0.5
  ) +
  theme_classic(
    base_size = 12
  ) +
  labs(
    title = "Hallmark pathway concordance",
    x = "NES: GSE197561",
    y = "NES: GSE111585"
  )

ggsave(
  filename = file.path(
    out_dir,
    "Part3_Hallmark_NES_Concordance_Scatter.pdf"
  ),
  plot = p_scatter,
  width = 7,
  height = 6.5
)

ggsave(
  filename = file.path(
    out_dir,
    "Part3_Hallmark_NES_Concordance_Scatter.png"
  ),
  plot = p_scatter,
  width = 7,
  height = 6.5,
  dpi = 300
)

heat_df <- concordance[
  is.finite(
    concordance$NES_GSE197561
  ) &
    is.finite(
      concordance$NES_GSE111585
    ),
  ,
  drop = FALSE
]

nes_mat <- heat_df[
  ,
  c(
    "NES_GSE197561",
    "NES_GSE111585"
  )
]

rownames(nes_mat) <-
  heat_df$Pathway

nes_mat <- as.matrix(
  nes_mat
)

colnames(nes_mat) <- c(
  "GSE197561",
  "GSE111585"
)

pdf(
  file.path(
    out_dir,
    "Part3_Hallmark_NES_Heatmap_All.pdf"
  ),
  width = 6,
  height = 11
)

pheatmap::pheatmap(
  nes_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Hallmark NES"
)

dev.off()

step2_results <- list(
  GMT = gmt_file,
  GSE197561 = gsea197,
  GSE111585 = gsea111,
  concordance = concordance,
  strict_concordant = strict_concordant,
  supported_concordant = supported_concordant,
  key_pathways = key_table
)

saveRDS(
  step2_results,
  file = file.path(
    out_dir,
    "Part3_Step2_Hallmark_V2_COMPLETE.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step2_Hallmark_V2_sessionInfo.txt"
  )
)

cat("\n\n")
cat("####################################################\n")
cat("STEP 2 V2 HALLMARK GSEA FINISHED\n")
cat("####################################################\n\n")

cat(
  "GSE197561 significant Hallmarks =",
  nrow(gsea197$sig),
  "\n"
)

cat(
  "GSE111585 significant Hallmarks =",
  nrow(gsea111$sig),
  "\n"
)

cat(
  "Strict concordant FDR-both =",
  nrow(strict_concordant),
  "\n"
)

cat(
  "Supported concordant =",
  nrow(supported_concordant),
  "\n"
)

cat(
  "GSE197561 unique warnings =",
  length(gsea197$warnings),
  "\n"
)

cat(
  "GSE111585 unique warnings =",
  length(gsea111$warnings),
  "\n"
)

cat(
  "\nResults saved to:\n",
  out_dir,
  "\n"
)
