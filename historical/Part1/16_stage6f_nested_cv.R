rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(glmnet)
library(timeROC)
library(survminer)
library(dplyr)

set.seed(20260813)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

gene_module <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.rds"
)

deg_df <- readRDS(
  "TCGA_OSCC_Nature8_DEG_results.rds"
)

clinical_wgcna <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_clinical_aligned.rds"
)

deg_ids <- unique(
  deg_df$gene_id[
    deg_df$DEG_status %in% c("Up", "Down")
  ]
)

candidate_table <- gene_module %>%
  dplyr::filter(
    module_color != "grey",
    gene_id %in% deg_ids
  ) %>%
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  )

deg_small <- deg_df %>%
  dplyr::select(
    gene_id,
    log2FoldChange,
    padj,
    DEG_status
  )

candidate_table <- dplyr::left_join(
  candidate_table,
  deg_small,
  by = "gene_id"
)

candidate_ids <- candidate_table$gene_id

cat("\n========================================\n")
cat("Candidate universe\n")
cat("========================================\n")

cat(
  "DEGs:",
  length(deg_ids),
  "\n"
)

cat(
  "DEG × non-grey WGCNA genes:",
  length(candidate_ids),
  "\n"
)

cat("\nCandidate modules:\n")

print(
  sort(
    table(candidate_table$module_color),
    decreasing = TRUE
  )
)

surv_keep <- (
  !is.na(clinical_wgcna$OS_time) &
  clinical_wgcna$OS_time > 0 &
  !is.na(clinical_wgcna$OS_event)
)

clinical_surv <- clinical_wgcna[
  surv_keep,
  ,
  drop = FALSE
]

x_all <- as.matrix(
  datExpr[
    surv_keep,
    candidate_ids,
    drop = FALSE
  ]
)

n <- nrow(clinical_surv)

cat("\nPatients =", n, "\n")
cat("Events:\n")

print(
  table(clinical_surv$OS_event)
)

outer_folds   <- 10
outer_repeats <- 10
inner_folds   <- 10

total_models <- outer_folds * outer_repeats

oof_list <- vector(
  "list",
  total_models
)

model_list <- vector(
  "list",
  total_models
)

selection_list <- vector(
  "list",
  total_models
)

counter <- 1

screen_unicox <- function(
    x_train,
    time,
    event,
    p_cutoff = 0.05) {

  genes <- colnames(x_train)

  pvals <- rep(
    NA_real_,
    length(genes)
  )

  names(pvals) <- genes

  for (j in seq_along(genes)) {

    z <- x_train[, j]

    fit <- try(
      survival::coxph(
        survival::Surv(
          time,
          event
        ) ~ z
      ),
      silent = TRUE
    )

    if (!inherits(fit, "try-error")) {

      sm <- summary(fit)

      pvals[j] <-
        sm$coefficients[
          1,
          "Pr(>|z|)"
        ]
    }
  }

  selected <- names(pvals)[
    !is.na(pvals) &
      pvals < p_cutoff
  ]

  list(
    selected = selected,
    pvalues = pvals
  )
}

for (rep_id in seq_len(outer_repeats)) {

  set.seed(
    20260813 + rep_id
  )

  outer_id <- rep(
    NA_integer_,
    n
  )

  alive_idx <- which(
    clinical_surv$OS_event == 0
  )

  dead_idx <- which(
    clinical_surv$OS_event == 1
  )

  outer_id[alive_idx] <- sample(
    rep(
      seq_len(outer_folds),
      length.out = length(alive_idx)
    )
  )

  outer_id[dead_idx] <- sample(
    rep(
      seq_len(outer_folds),
      length.out = length(dead_idx)
    )
  )

  for (fold_id in seq_len(outer_folds)) {

    test_idx <- which(
      outer_id == fold_id
    )

    train_idx <- which(
      outer_id != fold_id
    )

    x_train_all <- x_all[
      train_idx,
      ,
      drop = FALSE
    ]

    x_test_all <- x_all[
      test_idx,
      ,
      drop = FALSE
    ]

    time_train <-
      clinical_surv$OS_time[
        train_idx
      ]

    event_train <-
      clinical_surv$OS_event[
        train_idx
      ]

    uni <- screen_unicox(
      x_train = x_train_all,
      time = time_train,
      event = event_train,
      p_cutoff = 0.05
    )

    uni_genes <- uni$selected
    uni_n <- length(uni_genes)

    final_genes <- character(0)
    final_beta <- numeric(0)

    lambda_min <- NA_real_

    if (uni_n >= 2) {

      x_train <- x_train_all[
        ,
        uni_genes,
        drop = FALSE
      ]

      x_test <- x_test_all[
        ,
        uni_genes,
        drop = FALSE
      ]

      inner_id <- rep(
        NA_integer_,
        length(train_idx)
      )

      inner_alive <- which(
        event_train == 0
      )

      inner_dead <- which(
        event_train == 1
      )

      inner_id[inner_alive] <- sample(
        rep(
          seq_len(inner_folds),
          length.out = length(inner_alive)
        )
      )

      inner_id[inner_dead] <- sample(
        rep(
          seq_len(inner_folds),
          length.out = length(inner_dead)
        )
      )

      y_train <- survival::Surv(
        time_train,
        event_train
      )

      cvfit <- try(
        glmnet::cv.glmnet(
          x = x_train,
          y = y_train,
          family = "cox",
          alpha = 1,
          foldid = inner_id,
          standardize = TRUE,
          type.measure = "deviance"
        ),
        silent = TRUE
      )

      if (!inherits(
        cvfit,
        "try-error"
      )) {

        lambda_min <-
          cvfit$lambda.min

        beta <- as.numeric(
          coef(
            cvfit,
            s = "lambda.min"
          )
        )

        names(beta) <-
          colnames(x_train)

        final_genes <- names(beta)[
          beta != 0
        ]

        final_beta <- beta[
          final_genes
        ]
      }
    }

    if (uni_n == 1) {

      gene <- uni_genes[1]

      tmp_train <- data.frame(
        time = time_train,
        event = event_train,
        expression =
          x_train_all[
            ,
            gene
          ]
      )

      one_fit <- try(
        survival::coxph(
          survival::Surv(
            time,
            event
          ) ~ expression,
          data = tmp_train
        ),
        silent = TRUE
      )

      if (!inherits(
        one_fit,
        "try-error"
      )) {

        final_genes <- gene

        final_beta <-
          stats::coef(
            one_fit
          )[1]

        names(final_beta) <-
          gene
      }
    }

    final_n <- length(
      final_genes
    )

    if (final_n == 0) {

      risk_test_raw <- rep(
        0,
        length(test_idx)
      )

      risk_test_z <- rep(
        0,
        length(test_idx)
      )

    } else {

      lp_train <- as.numeric(
        x_train_all[
          ,
          final_genes,
          drop = FALSE
        ] %*%
          final_beta[
            final_genes
          ]
      )

      risk_test_raw <- as.numeric(
        x_test_all[
          ,
          final_genes,
          drop = FALSE
        ] %*%
          final_beta[
            final_genes
          ]
      )

      train_lp_mean <- mean(
        lp_train
      )

      train_lp_sd <- sd(
        lp_train
      )

      if (
        is.na(train_lp_sd) ||
        train_lp_sd == 0
      ) {

        risk_test_z <- rep(
          0,
          length(test_idx)
        )

      } else {

        risk_test_z <- (
          risk_test_raw -
            train_lp_mean
        ) /
          train_lp_sd
      }
    }

    oof_list[[counter]] <- data.frame(
      repeat_id = rep_id,
      fold_id = fold_id,

      patient_barcode =
        clinical_surv$patient_barcode[
          test_idx
        ],

      OS_time =
        clinical_surv$OS_time[
          test_idx
        ],

      OS_event =
        clinical_surv$OS_event[
          test_idx
        ],

      risk_raw =
        risk_test_raw,

      risk_score =
        risk_test_z,

      uni_gene_n =
        uni_n,

      final_gene_n =
        final_n,

      stringsAsFactors = FALSE
    )

    model_list[[counter]] <- data.frame(
      repeat_id = rep_id,
      fold_id = fold_id,

      uni_gene_n =
        uni_n,

      final_gene_n =
        final_n,

      lambda_min =
        lambda_min,

      null_model =
        final_n == 0,

      stringsAsFactors = FALSE
    )

    if (final_n > 0) {

      selection_list[[counter]] <- data.frame(
        repeat_id = rep_id,
        fold_id = fold_id,
        gene_id = final_genes,
        coefficient =
          as.numeric(
            final_beta[
              final_genes
            ]
          ),
        stringsAsFactors = FALSE
      )

    } else {

      selection_list[[counter]] <-
        data.frame(
          repeat_id = integer(0),
          fold_id = integer(0),
          gene_id = character(0),
          coefficient = numeric(0)
        )
    }

    counter <- counter + 1
  }

  cat(
    "Completed repeat",
    rep_id,
    "/",
    outer_repeats,
    "\n"
  )
}

oof_all <- dplyr::bind_rows(
  oof_list
)

model_summary <- dplyr::bind_rows(
  model_list
)

selection_all <- dplyr::bind_rows(
  selection_list
)

oof_patient <- oof_all %>%
  dplyr::group_by(
    patient_barcode
  ) %>%
  dplyr::summarise(

    OS_time =
      dplyr::first(OS_time),

    OS_event =
      dplyr::first(OS_event),

    risk_score =
      mean(
        risk_score,
        na.rm = TRUE
      ),

    risk_sd =
      sd(
        risk_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

stopifnot(
  nrow(oof_patient) ==
    nrow(clinical_surv)
)

c_oof <- survival::concordance(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_score,
  data = oof_patient,
  reverse = TRUE
)

roc_oof <- timeROC::timeROC(
  T = oof_patient$OS_time,
  delta = oof_patient$OS_event,
  marker = oof_patient$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(
    365,
    1095,
    1825
  ),
  iid = TRUE
)

risk_cutoff <- median(
  oof_patient$risk_score,
  na.rm = TRUE
)

oof_patient$risk_group <- ifelse(
  oof_patient$risk_score >=
    risk_cutoff,
  "High",
  "Low"
)

oof_patient$risk_group <- factor(
  oof_patient$risk_group,
  levels = c(
    "Low",
    "High"
  )
)

fit_km <- survival::survfit(
  survival::Surv(
    OS_time,
    OS_event
  ) ~ risk_group,
  data = oof_patient
)

km_plot <- survminer::ggsurvplot(
  fit_km,
  data = oof_patient,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title =
    "Nested-CV out-of-fold prognostic prediction",
  legend.title =
    "Risk group"
)

pdf(
  "Figure_Nature8_OSCC_Stage6F_NestedCV_KM.pdf",
  width = 7,
  height = 7
)

print(km_plot)

dev.off()

pdf(
  "Figure_Nature8_OSCC_Stage6F_NestedCV_timeROC.pdf",
  width = 7,
  height = 6
)

plot(
  roc_oof,
  time = 365,
  lwd = 2,
  title = FALSE
)

plot(
  roc_oof,
  time = 1095,
  add = TRUE,
  lwd = 2
)

plot(
  roc_oof,
  time = 1825,
  add = TRUE,
  lwd = 2
)

legend(
  "bottomright",
  legend = c(
    paste0(
      "1-year AUC = ",
      round(
        roc_oof$AUC[1],
        3
      )
    ),

    paste0(
      "3-year AUC = ",
      round(
        roc_oof$AUC[2],
        3
      )
    ),

    paste0(
      "5-year AUC = ",
      round(
        roc_oof$AUC[3],
        3
      )
    )
  ),
  lwd = 2
)

title(
  "Nested-CV out-of-fold ROC"
)

dev.off()

if (
  nrow(selection_all) > 0
) {

  nested_selection <- selection_all %>%
    dplyr::count(
      gene_id,
      name = "selection_count"
    )

  nested_selection$selection_frequency <-
    nested_selection$selection_count /
    total_models

  nested_selection <- dplyr::left_join(
    nested_selection,
    candidate_table %>%
      dplyr::select(
        gene_id,
        gene_name,
        module_color,
        log2FoldChange,
        DEG_status
      ),
    by = "gene_id"
  )

  nested_selection <- nested_selection %>%
    dplyr::arrange(
      dplyr::desc(
        selection_frequency
      )
    )

} else {

  nested_selection <- data.frame()
}

write.table(
  candidate_table,
  "TCGA_OSCC_Nature8_Stage6F_DEG_nonGrey_candidates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  oof_all,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_AllPredictions.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  oof_patient,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_PatientPredictions.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  model_summary,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_ModelSummary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  nested_selection,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_GeneSelection.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  oof_patient,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_PatientPredictions.rds"
)

saveRDS(
  model_summary,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_ModelSummary.rds"
)

saveRDS(
  nested_selection,
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_GeneSelection.rds"
)

cat("\n========================================\n")
cat("STAGE 6F SUMMARY\n")
cat("========================================\n")

cat(
  "Candidate genes:",
  length(candidate_ids),
  "\n"
)

cat(
  "Outer models:",
  nrow(model_summary),
  "\n"
)

cat(
  "Null models:",
  sum(model_summary$null_model),
  "\n"
)

cat(
  "Null model percent:",
  round(
    mean(model_summary$null_model) *
      100,
    1
  ),
  "%\n"
)

cat(
  "Median Uni-Cox genes:",
  median(
    model_summary$uni_gene_n
  ),
  "\n"
)

cat(
  "Median final genes:",
  median(
    model_summary$final_gene_n
  ),
  "\n"
)

cat(
  "OOF C-index:",
  c_oof$concordance,
  "\n"
)

cat(
  "1-year AUC:",
  roc_oof$AUC[1],
  "\n"
)

cat(
  "3-year AUC:",
  roc_oof$AUC[2],
  "\n"
)

cat(
  "5-year AUC:",
  roc_oof$AUC[3],
  "\n"
)

cat("========================================\n")
