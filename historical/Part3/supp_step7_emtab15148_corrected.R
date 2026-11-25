rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required.")
}

suppressPackageStartupMessages({
  library(Matrix)
})

project_dir <- file.path(project_root, "Part3_Cisplatin")

sc_dir <- file.path(
  project_dir,
  "E-MTAB-15148"
)

old_step7_dir <- file.path(
  project_dir,
  "Step7_EMTAB15148_TemporalValidation"
)

preflight_file <- file.path(
  sc_dir,
  "E-MTAB-15148_Step7_Preflight.rds"
)

step3_file <- file.path(
  project_dir,
  "Step3_Gene_Convergence",
  "Part3_Step3_GeneConvergence_COMPLETE.rds"
)

old_step7_file <- file.path(
  old_step7_dir,
  "Part3_Step7_EMTAB15148_COMPLETE.rds"
)

required_files <- c(
  preflight_file,
  step3_file,
  old_step7_file
)

if (!all(file.exists(required_files))) {
  cat("\nMissing files:\n")
  print(required_files[!file.exists(required_files)])
  stop("Required file missing.")
}

out_dir <- file.path(
  project_dir,
  "Step7_EMTAB15148_TemporalValidation_CORRECTED"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

preflight <- readRDS(
  preflight_file
)

step3 <- readRDS(
  step3_file
)

old_step7 <- readRDS(
  old_step7_file
)

expr <- preflight$expr

features <- as.data.frame(
  preflight$features
)

barcodes <- as.character(
  preflight$barcodes
)

md <- as.data.frame(
  old_step7$demux_metadata
)

required_md_cols <- c(
  "Barcode",
  "HTO_Global",
  "HTO_MaxID"
)

if (!all(required_md_cols %in% colnames(md))) {
  stop(
    "Saved Step 7 metadata does not contain required HTO columns."
  )
}

cat("\n========== ACTUAL HTO LABELS ==========\n")

print(
  table(
    md$HTO_MaxID,
    md$HTO_Global,
    useNA = "ifany"
  )
)

state_map <- c(
  "PDC1" = "PDC1_excluded",
  "SCC9" = "CisS",
  "CISR-4M" = "CisR4M",
  "CISR-8M" = "CisR8M"
)

md$State <- unname(
  state_map[
    as.character(md$HTO_MaxID)
  ]
)

singlet_tags <- unique(
  as.character(
    md$HTO_MaxID[
      md$HTO_Global == "Singlet"
    ]
  )
)

unmapped_singlet_tags <- setdiff(
  singlet_tags,
  names(state_map)
)

if (length(unmapped_singlet_tags) > 0) {
  stop(
    paste0(
      "Unmapped Singlet HTO labels: ",
      paste(
        unmapped_singlet_tags,
        collapse = ", "
      )
    )
  )
}

keep_states <- c(
  "CisS",
  "CisR4M",
  "CisR8M"
)

keep <- md$HTO_Global == "Singlet" &
  md$State %in% keep_states

final_cells <- md$Barcode[
  keep
]

state_vector <- md$State[
  keep
]

names(state_vector) <- final_cells

state_vector <- factor(
  state_vector,
  levels = keep_states
)

state_counts <- data.frame(
  State = keep_states,
  N_cells = as.integer(
    table(state_vector)[keep_states]
  ),
  stringsAsFactors = FALSE
)

cat("\n========== CORRECTED STATE COUNTS ==========\n")

print(state_counts)

expected_counts <- c(
  CisS = 1645,
  CisR4M = 1073,
  CisR8M = 1743
)

observed_counts <- setNames(
  state_counts$N_cells,
  state_counts$State
)

if (!all(
  observed_counts[names(expected_counts)] ==
    expected_counts
)) {
  stop(
    "Corrected state counts do not match the existing HTODemux Singlet counts."
  )
}

cat(
  "\nTotal SCC9 singlets =",
  length(final_cells),
  "\n"
)

feature_name <- as.character(
  features[, 2]
)

feature_type <- as.character(
  features[, 3]
)

gex_idx <- which(
  feature_type == "Gene Expression"
)

gex_counts <- expr[
  gex_idx,
  ,
  drop = FALSE
]

gex_symbols_raw <- feature_name[
  gex_idx
]

colnames(gex_counts) <- barcodes

lib_size <- Matrix::colSums(
  gex_counts[
    ,
    final_cells,
    drop = FALSE
  ]
)

detected_genes <- Matrix::colSums(
  gex_counts[
    ,
    final_cells,
    drop = FALSE
  ] > 0
)

rna_qc <- data.frame(
  Barcode = final_cells,
  State = as.character(
    state_vector[
      final_cells
    ]
  ),
  UMI = as.numeric(
    lib_size[
      final_cells
    ]
  ),
  Detected_Genes = as.numeric(
    detected_genes[
      final_cells
    ]
  ),
  stringsAsFactors = FALSE
)

priority_table <- as.data.frame(
  step3$priority_pool
)

ultra_table <- as.data.frame(
  step3$ultra_strict
)

priority_genes <- unique(
  as.character(
    priority_table$Gene
  )
)

ultra_genes <- unique(
  as.character(
    ultra_table$Gene
  )
)

target_genes <- unique(
  c(
    priority_genes,
    ultra_genes,
    "BASP1"
  )
)

collapse_target_counts <- function(
    sparse_matrix,
    raw_symbols,
    genes,
    cells
) {

  out <- matrix(
    0,
    nrow = length(genes),
    ncol = length(cells),
    dimnames = list(
      genes,
      cells
    )
  )

  for (i in seq_along(genes)) {

    gene <- genes[i]

    idx <- which(
      raw_symbols == gene
    )

    if (length(idx) > 0) {

      out[i, ] <- Matrix::colSums(
        sparse_matrix[
          idx,
          cells,
          drop = FALSE
        ]
      )
    }
  }

  out
}

target_counts <- collapse_target_counts(
  sparse_matrix = gex_counts,
  raw_symbols = gex_symbols_raw,
  genes = target_genes,
  cells = final_cells
)

genes_detected_any <- rownames(
  target_counts
)[
  rowSums(target_counts) > 0
]

priority_available <- intersect(
  priority_genes,
  genes_detected_any
)

ultra_available <- intersect(
  ultra_genes,
  genes_detected_any
)

cat(
  "\nPriority available =",
  length(priority_available),
  "/",
  length(priority_genes),
  "\n"
)

cat(
  "Ultra available =",
  length(ultra_available),
  "/",
  length(ultra_genes),
  "\n"
)

target_norm <- sweep(
  target_counts,
  2,
  as.numeric(
    lib_size[
      colnames(target_counts)
    ]
  ),
  FUN = "/"
)

target_norm <- target_norm * 10000

target_log <- log1p(
  target_norm
)

make_signed_score <- function(
    log_expression,
    gene_table,
    genes
) {

  genes <- intersect(
    genes,
    rownames(log_expression)
  )

  genes <- genes[
    rowSums(
      log_expression[
        genes,
        ,
        drop = FALSE
      ]
    ) > 0
  ]

  x <- log_expression[
    genes,
    ,
    drop = FALSE
  ]

  gene_sd <- apply(
    x,
    1,
    sd,
    na.rm = TRUE
  )

  good <- is.finite(gene_sd) &
    gene_sd > 0

  x <- x[
    good,
    ,
    drop = FALSE
  ]

  if (nrow(x) < 2) {
    stop(
      "Too few genes available for signature score."
    )
  }

  x_z <- t(
    scale(
      t(x)
    )
  )

  meta_z <- gene_table$Meta_Z[
    match(
      rownames(x_z),
      gene_table$Gene
    )
  ]

  direction_sign <- ifelse(
    meta_z >= 0,
    1,
    -1
  )

  signed_x <- sweep(
    x_z,
    1,
    direction_sign,
    FUN = "*"
  )

  score <- colMeans(
    signed_x,
    na.rm = TRUE
  )

  names(score) <- colnames(
    signed_x
  )

  list(
    score = score,
    genes_used = rownames(signed_x)
  )
}

priority_score_result <- make_signed_score(
  log_expression = target_log,
  gene_table = priority_table,
  genes = priority_genes
)

ultra_score_result <- make_signed_score(
  log_expression = target_log,
  gene_table = ultra_table,
  genes = ultra_genes
)

cell_scores <- data.frame(
  Barcode = final_cells,
  State = as.character(
    state_vector[
      final_cells
    ]
  ),
  Priority_Score =
    priority_score_result$score[
      final_cells
    ],
  Ultra_Score =
    ultra_score_result$score[
      final_cells
    ],
  stringsAsFactors = FALSE
)

summarize_score <- function(
    score_table,
    score_column
) {

  output <- data.frame()

  for (state in keep_states) {

    x <- score_table[
      score_table$State == state,
      score_column
    ]

    output <- rbind(
      output,
      data.frame(
        Score = score_column,
        State = state,
        N_cells = length(x),
        Mean = mean(
          x,
          na.rm = TRUE
        ),
        Median = median(
          x,
          na.rm = TRUE
        ),
        SD = sd(
          x,
          na.rm = TRUE
        ),
        Q25 = unname(
          quantile(
            x,
            0.25,
            na.rm = TRUE
          )
        ),
        Q75 = unname(
          quantile(
            x,
            0.75,
            na.rm = TRUE
          )
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  output
}

score_summary <- rbind(
  summarize_score(
    cell_scores,
    "Priority_Score"
  ),
  summarize_score(
    cell_scores,
    "Ultra_Score"
  )
)

score_trajectory <- data.frame()

for (score_name in c(
  "Priority_Score",
  "Ultra_Score"
)) {

  x <- score_summary[
    score_summary$Score == score_name,
    ,
    drop = FALSE
  ]

  mean_s <- x$Mean[
    x$State == "CisS"
  ]

  mean_4 <- x$Mean[
    x$State == "CisR4M"
  ]

  mean_8 <- x$Mean[
    x$State == "CisR8M"
  ]

  median_s <- x$Median[
    x$State == "CisS"
  ]

  median_4 <- x$Median[
    x$State == "CisR4M"
  ]

  median_8 <- x$Median[
    x$State == "CisR8M"
  ]

  score_trajectory <- rbind(
    score_trajectory,
    data.frame(
      Score = score_name,

      CisS_Mean = mean_s,
      CisR4M_Mean = mean_4,
      CisR8M_Mean = mean_8,

      Early_minus_Sensitive =
        mean_4 - mean_s,

      Late_minus_Sensitive =
        mean_8 - mean_s,

      Late_minus_Early =
        mean_8 - mean_4,

      Mean_Monotonic_Increase =
        mean_s <= mean_4 &&
        mean_4 <= mean_8,

      CisS_Median = median_s,
      CisR4M_Median = median_4,
      CisR8M_Median = median_8,

      Median_Monotonic_Increase =
        median_s <= median_4 &&
        median_4 <= median_8,

      stringsAsFactors = FALSE
    )
  )
}

gene_trajectory_summary <- function(
    genes,
    gene_table,
    log_expression,
    raw_counts,
    states
) {

  output <- data.frame()

  for (gene in genes) {

    if (!(gene %in% rownames(log_expression))) {
      next
    }

    x <- log_expression[
      gene,
      ,
      drop = TRUE
    ]

    raw_x <- raw_counts[
      gene,
      ,
      drop = TRUE
    ]

    mean_s <- mean(
      x[
        states == "CisS"
      ],
      na.rm = TRUE
    )

    mean_4 <- mean(
      x[
        states == "CisR4M"
      ],
      na.rm = TRUE
    )

    mean_8 <- mean(
      x[
        states == "CisR8M"
      ],
      na.rm = TRUE
    )

    det_s <- mean(
      raw_x[
        states == "CisS"
      ] > 0
    )

    det_4 <- mean(
      raw_x[
        states == "CisR4M"
      ] > 0
    )

    det_8 <- mean(
      raw_x[
        states == "CisR8M"
      ] > 0
    )

    meta_z <- gene_table$Meta_Z[
      match(
        gene,
        gene_table$Gene
      )
    ]

    expected_sign <- ifelse(
      meta_z >= 0,
      1,
      -1
    )

    signed_s <- expected_sign * mean_s
    signed_4 <- expected_sign * mean_4
    signed_8 <- expected_sign * mean_8

    output <- rbind(
      output,
      data.frame(
        Gene = gene,
        Meta_Z = meta_z,

        Expected_Direction =
          ifelse(
            expected_sign > 0,
            "Higher_in_resistance",
            "Lower_in_resistance"
          ),

        CisS_Mean = mean_s,
        CisR4M_Mean = mean_4,
        CisR8M_Mean = mean_8,

        Early_minus_Sensitive =
          mean_4 - mean_s,

        Late_minus_Sensitive =
          mean_8 - mean_s,

        Late_minus_Early =
          mean_8 - mean_4,

        Signed_Early_vs_Sensitive =
          expected_sign *
          (mean_4 - mean_s),

        Signed_Late_vs_Sensitive =
          expected_sign *
          (mean_8 - mean_s),

        Supports_Early =
          signed_4 > signed_s,

        Supports_Late =
          signed_8 > signed_s,

        Monotonic_Expected_Direction =
          signed_s <= signed_4 &&
          signed_4 <= signed_8,

        Detection_CisS = det_s,
        Detection_CisR4M = det_4,
        Detection_CisR8M = det_8,

        stringsAsFactors = FALSE
      )
    )
  }

  output
}

states_for_matrix <- as.character(
  state_vector[
    colnames(target_log)
  ]
)

priority_gene_trajectory <- gene_trajectory_summary(
  genes = priority_score_result$genes_used,
  gene_table = priority_table,
  log_expression = target_log,
  raw_counts = target_counts,
  states = states_for_matrix
)

ultra_gene_trajectory <- gene_trajectory_summary(
  genes = ultra_score_result$genes_used,
  gene_table = ultra_table,
  log_expression = target_log,
  raw_counts = target_counts,
  states = states_for_matrix
)

ultra_gene_trajectory <- ultra_gene_trajectory[
  order(
    -as.integer(
      ultra_gene_trajectory$Monotonic_Expected_Direction
    ),
    -ultra_gene_trajectory$Signed_Late_vs_Sensitive
  ),
  ,
  drop = FALSE
]

basp1_trajectory <- data.frame()

if ("BASP1" %in% rownames(target_log)) {

  x <- target_log[
    "BASP1",
    ,
    drop = TRUE
  ]

  raw_x <- target_counts[
    "BASP1",
    ,
    drop = TRUE
  ]

  basp1_trajectory <- data.frame(
    State = keep_states,

    Mean_logNorm = c(
      mean(
        x[
          states_for_matrix ==
            "CisS"
        ],
        na.rm = TRUE
      ),

      mean(
        x[
          states_for_matrix ==
            "CisR4M"
        ],
        na.rm = TRUE
      ),

      mean(
        x[
          states_for_matrix ==
            "CisR8M"
        ],
        na.rm = TRUE
      )
    ),

    Detection_Fraction = c(
      mean(
        raw_x[
          states_for_matrix ==
            "CisS"
        ] > 0
      ),

      mean(
        raw_x[
          states_for_matrix ==
            "CisR4M"
        ] > 0
      ),

      mean(
        raw_x[
          states_for_matrix ==
            "CisR8M"
        ] > 0
      )
    ),

    stringsAsFactors = FALSE
  )
}

gene_direction_summary <- data.frame(
  Gene_Set = c(
    "Priority",
    "Ultra"
  ),

  N_genes = c(
    nrow(priority_gene_trajectory),
    nrow(ultra_gene_trajectory)
  ),

  N_support_early = c(
    sum(
      priority_gene_trajectory$Supports_Early,
      na.rm = TRUE
    ),
    sum(
      ultra_gene_trajectory$Supports_Early,
      na.rm = TRUE
    )
  ),

  N_support_late = c(
    sum(
      priority_gene_trajectory$Supports_Late,
      na.rm = TRUE
    ),
    sum(
      ultra_gene_trajectory$Supports_Late,
      na.rm = TRUE
    )
  ),

  N_monotonic = c(
    sum(
      priority_gene_trajectory$Monotonic_Expected_Direction,
      na.rm = TRUE
    ),
    sum(
      ultra_gene_trajectory$Monotonic_Expected_Direction,
      na.rm = TRUE
    )
  ),

  stringsAsFactors = FALSE
)

write.table(
  md,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_HTO_Metadata.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  state_counts,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_StateCounts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  rna_qc,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_RNA_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cell_scores,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_CellScores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  score_summary,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_ScoreSummary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  score_trajectory,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_ScoreTrajectory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  priority_gene_trajectory,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_Priority52_Trajectory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ultra_gene_trajectory,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_Ultra7_Trajectory.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_direction_summary,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_GeneDirectionSummary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  basp1_trajectory,
  file = file.path(
    out_dir,
    "EMTAB15148_CORRECTED_BASP1.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

step7_corrected <- list(
  demux_metadata = md,
  state_counts = state_counts,
  rna_qc = rna_qc,
  priority_genes_used =
    priority_score_result$genes_used,
  ultra_genes_used =
    ultra_score_result$genes_used,
  cell_scores = cell_scores,
  score_summary = score_summary,
  score_trajectory = score_trajectory,
  priority_gene_trajectory =
    priority_gene_trajectory,
  ultra_gene_trajectory =
    ultra_gene_trajectory,
  gene_direction_summary =
    gene_direction_summary,
  basp1_trajectory =
    basp1_trajectory
)

saveRDS(
  step7_corrected,
  file = file.path(
    out_dir,
    "Part3_Step7_EMTAB15148_CORRECTED_COMPLETE.rds"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    out_dir,
    "Part3_Step7_CORRECTED_sessionInfo.txt"
  )
)

cat("\n")
cat("####################################################\n")
cat("STEP 7 CORRECTION FINISHED SUCCESSFULLY\n")
cat("####################################################\n\n")

print(state_counts)
print(score_trajectory)
print(gene_direction_summary)
print(ultra_gene_trajectory)
print(basp1_trajectory)
