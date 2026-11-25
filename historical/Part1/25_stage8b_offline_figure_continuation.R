project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)
library(survival)
library(timeROC)
library(grid)

setwd(project_root)

gsea_top <- data.frame(

  Description = c(
    "E2F TARGETS",
    "G2M CHECKPOINT",
    "EPITHELIAL MESENCHYMAL TRANSITION",
    "OXIDATIVE PHOSPHORYLATION",
    "INTERFERON GAMMA RESPONSE",
    "INTERFERON ALPHA RESPONSE",
    "MYOGENESIS",
    "ADIPOGENESIS",
    "FATTY ACID METABOLISM",
    "MITOTIC SPINDLE"
  ),

  NES = c(
     3.497424,
     3.278074,
     2.874844,
    -2.665429,
     2.713764,
     3.042602,
    -2.503257,
    -2.482261,
    -2.425878,
     2.281158
  ),

  padj = c(
    3.333533e-44,
    1.515064e-34,
    1.098413e-22,
    6.906693e-19,
    8.625290e-19,
    9.510290e-19,
    3.491964e-16,
    8.255908e-16,
    2.246174e-12,
    2.416417e-11
  )
)

gsea_top$Direction <- ifelse(
  gsea_top$NES > 0,
  "Tumor enriched",
  "Normal enriched"
)

gsea_top <- gsea_top %>%
  dplyr::arrange(NES)

gsea_top$Description <- factor(
  gsea_top$Description,
  levels = gsea_top$Description
)

p_gsea <- ggplot(
  gsea_top,
  aes(
    x = NES,
    y = Description
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = NES,
      yend = Description,
      linetype = Direction
    ),
    linewidth = 0.9
  ) +

  geom_point(
    aes(
      shape = Direction
    ),
    size = 3.2
  ) +

  labs(
    title = "Hallmark pathway enrichment in OSCC",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    shape = NULL,
    linetype = NULL
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    axis.text.y = element_text(
      size = 10
    ),

    legend.position = "bottom"
  )

ggsave(
  "Figure1D_OSCC_Hallmark_GSEA.pdf",
  p_gsea,
  width = 9.5,
  height = 6.5
)

ggsave(
  "Figure1D_OSCC_Hallmark_GSEA.png",
  p_gsea,
  width = 9.5,
  height = 6.5,
  dpi = 600
)

module_cox <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_Module_OS_Cox_COLOR.rds"
)

module_cox$module <- gsub(
  "^ME",
  "",
  module_cox$Module
)

module_cox <- module_cox %>%
  dplyr::arrange(HR)

module_cox$module <- factor(
  module_cox$module,
  levels = module_cox$module
)

module_cox$label <- paste0(
  sprintf("%.2f", module_cox$HR),
  " (",
  sprintf("%.2f", module_cox$lower95),
  "–",
  sprintf("%.2f", module_cox$upper95),
  ")"
)

p_module_cox <- ggplot(
  module_cox,
  aes(
    x = HR,
    y = module
  )
) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +

  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0.16,
    linewidth = 0.7
  ) +

  geom_point(
    size = 3
  ) +

  scale_x_log10() +

  labs(
    title = "Module-level association with overall survival",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "Figure2D_WGCNA_Module_OS_Cox.pdf",
  p_module_cox,
  width = 8,
  height = 5.8
)

ggsave(
  "Figure2D_WGCNA_Module_OS_Cox.png",
  p_module_cox,
  width = 8,
  height = 5.8,
  dpi = 600
)

nested_selection <- readRDS(
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_GeneSelection.rds"
)

gene_module <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.rds"
)

deg_df <- readRDS(
  "TCGA_OSCC_Nature8_DEG_results.rds"
)

n_deg <- nrow(
  deg_df
)

n_non_grey <- sum(
  gene_module$module_color != "grey"
)

candidate_non_grey_deg <- gene_module %>%
  dplyr::filter(
    module_color != "grey",
    gene_id %in% deg_df$gene_id
  ) %>%
  dplyr::distinct(
    gene_id
  )

n_candidate <- nrow(
  candidate_non_grey_deg
)

n_stable8 <- sum(
  nested_selection$selection_frequency >= 0.80
)

flow_df <- data.frame(

  x = 1:5,

  label = c(
    paste0(
      "Differentially expressed genes\nn = ",
      n_deg
    ),

    paste0(
      "Non-grey WGCNA genes\nn = ",
      n_non_grey
    ),

    paste0(
      "DEG ∩ WGCNA\nn = ",
      n_candidate
    ),

    "Training-fold Cox + LASSO\nRepeated nested CV",

    paste0(
      "Stable 8-gene signature\nn = ",
      n_stable8
    )
  )
)

p_flow <- ggplot() +

  geom_segment(
    data = flow_df[-nrow(flow_df), ],
    aes(
      x = x + 0.31,
      xend = x + 0.69,
      y = 1,
      yend = 1
    ),
    arrow = arrow(
      length = unit(
        0.18,
        "cm"
      )
    ),
    linewidth = 0.8
  ) +

  geom_label(
    data = flow_df,
    aes(
      x = x,
      y = 1,
      label = label
    ),
    size = 3.4,
    label.padding = unit(
      0.38,
      "lines"
    )
  ) +

  coord_cartesian(
    xlim = c(
      0.5,
      5.5
    ),
    ylim = c(
      0.72,
      1.28
    ),
    clip = "off"
  ) +

  labs(
    title = "Prognostic feature-selection strategy"
  ) +

  theme_void() +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    )
  )

ggsave(
  "Figure3A_Prognostic_FeatureSelection_Flow.pdf",
  p_flow,
  width = 14,
  height = 3.8
)

ggsave(
  "Figure3A_Prognostic_FeatureSelection_Flow.png",
  p_flow,
  width = 14,
  height = 3.8,
  dpi = 600
)

stable8 <- nested_selection %>%
  dplyr::filter(
    selection_frequency >= 0.80
  ) %>%
  dplyr::arrange(
    selection_frequency
  )

stable8$gene_name <- factor(
  stable8$gene_name,
  levels = stable8$gene_name
)

p_stable8 <- ggplot(
  stable8,
  aes(
    x = selection_frequency,
    y = gene_name
  )
) +

  geom_col(
    width = 0.65
  ) +

  geom_vline(
    xintercept = 0.80,
    linetype = "dashed",
    linewidth = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        selection_frequency
      )
    ),
    hjust = -0.2,
    size = 3.8
  ) +

  scale_x_continuous(
    limits = c(
      0,
      1.08
    ),
    breaks = seq(
      0,
      1,
      0.2
    )
  ) +

  labs(
    title = "Stable genes selected in nested cross-validation",
    subtitle = "Selection frequency ≥ 0.80",
    x = "Selection frequency",
    y = NULL
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    axis.text.y = element_text(
      face = "italic"
    )
  )

ggsave(
  "Figure3C_NestedCV_Final8_Selection.pdf",
  p_stable8,
  width = 8,
  height = 5.8
)

ggsave(
  "Figure3C_NestedCV_Final8_Selection.png",
  p_stable8,
  width = 8,
  height = 5.8,
  dpi = 600
)

tcga <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_RiskScore.rds"
)

tcga_plot <- tcga %>%
  dplyr::arrange(
    risk_score
  )

tcga_plot$Patient_rank <- seq_len(
  nrow(tcga_plot)
)

p_tcga_risk <- ggplot(
  tcga_plot,
  aes(
    x = Patient_rank,
    y = risk_score
  )
) +

  geom_point(
    aes(
      shape = risk_group
    ),
    size = 1.7
  ) +

  geom_hline(
    yintercept = median(
      tcga_plot$risk_score
    ),
    linetype = "dashed"
  ) +

  labs(
    title = "TCGA risk-score distribution",
    x = "Patients ordered by risk score",
    y = "Risk score",
    shape = "Risk group"
  ) +

  theme_classic(
    base_size = 12
  )

ggsave(
  "Figure4A_TCGA_RiskScore_Distribution.pdf",
  p_tcga_risk,
  width = 8,
  height = 5
)

p_tcga_surv <- ggplot(
  tcga_plot,
  aes(
    x = Patient_rank,
    y = OS_time
  )
) +

  geom_point(
    aes(
      shape = factor(
        OS_event
      )
    ),
    size = 1.8
  ) +

  labs(
    title = "TCGA survival-status distribution",
    x = "Patients ordered by risk score",
    y = "Overall survival time (days)",
    shape = "OS event"
  ) +

  theme_classic(
    base_size = 12
  )

ggsave(
  "Figure4B_TCGA_SurvivalStatus.pdf",
  p_tcga_surv,
  width = 8,
  height = 5
)

internal_perf <- data.frame(

  Evaluation = factor(
    c(
      "TCGA apparent",
      "Nested-CV OOF"
    ),
    levels = c(
      "TCGA apparent",
      "Nested-CV OOF"
    )
  ),

  C_index = c(
    0.6730039,
    0.5695233
  )
)

p_internal <- ggplot(
  internal_perf,
  aes(
    x = Evaluation,
    y = C_index
  )
) +

  geom_hline(
    yintercept = 0.5,
    linetype = "dotted"
  ) +

  geom_point(
    size = 4
  ) +

  geom_segment(
    aes(
      x = 1,
      xend = 2,
      y = C_index[1],
      yend = C_index[2]
    ),
    linewidth = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        C_index
      )
    ),
    vjust = -1,
    size = 4
  ) +

  coord_cartesian(
    ylim = c(
      0.50,
      0.72
    )
  ) +

  labs(
    title = "Internal validation of prognostic discrimination",
    x = NULL,
    y = "C-index"
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "Figure4E_TCGA_InternalValidation_Cindex.pdf",
  p_internal,
  width = 7,
  height = 5
)

geo <- readRDS(
  "GSE41613_Final8Gene_ExternalValidation.rds"
)

geo_plot <- geo %>%
  dplyr::arrange(
    risk_score
  )

geo_plot$Patient_rank <- seq_len(
  nrow(geo_plot)
)

p_geo_risk <- ggplot(
  geo_plot,
  aes(
    x = Patient_rank,
    y = risk_score
  )
) +

  geom_point(
    aes(
      shape = risk_group
    ),
    size = 2
  ) +

  labs(
    title = "GSE41613 risk-score distribution",
    x = "Patients ordered by risk score",
    y = "Risk score",
    shape = "Risk group"
  ) +

  theme_classic(
    base_size = 12
  )

ggsave(
  "Figure5A_GSE41613_RiskScore_Distribution.pdf",
  p_geo_risk,
  width = 8,
  height = 5
)

p_geo_surv <- ggplot(
  geo_plot,
  aes(
    x = Patient_rank,
    y = OS_time_days
  )
) +

  geom_point(
    aes(
      shape = factor(
        OS_event
      )
    ),
    size = 2
  ) +

  labs(
    title = "GSE41613 survival-status distribution",
    x = "Patients ordered by risk score",
    y = "Overall survival time (days)",
    shape = "OS event"
  ) +

  theme_classic(
    base_size = 12
  )

ggsave(
  "Figure5B_GSE41613_SurvivalStatus.pdf",
  p_geo_surv,
  width = 8,
  height = 5
)

fit_combined <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_CombinedClinicalModel.rds"
)

s_comb <- summary(
  fit_combined
)

clin_forest <- data.frame(

  Variable = rownames(
    s_comb$coefficients
  ),

  HR = s_comb$coefficients[
    ,
    "exp(coef)"
  ],

  lower95 = s_comb$conf.int[
    ,
    "lower .95"
  ],

  upper95 = s_comb$conf.int[
    ,
    "upper .95"
  ],

  P = s_comb$coefficients[
    ,
    "Pr(>|z|)"
  ],

  stringsAsFactors = FALSE
)

clin_forest$Label <- clin_forest$Variable

clin_forest$Label[
  clin_forest$Variable == "age"
] <- "Age"

clin_forest$Label[
  clin_forest$Variable == "sexmale"
] <- "Male sex"

clin_forest$Label[
  clin_forest$Variable == "stage4Stage II"
] <- "Stage II"

clin_forest$Label[
  clin_forest$Variable == "stage4Stage III"
] <- "Stage III"

clin_forest$Label[
  clin_forest$Variable == "stage4Stage IVA"
] <- "Stage IVA"

clin_forest$Label[
  clin_forest$Variable == "stage4Stage IVB"
] <- "Stage IVB"

clin_forest$Label[
  clin_forest$Variable == "risk_z"
] <- "8-gene risk score (per SD)"

clin_forest$Label <- factor(
  clin_forest$Label,
  levels = rev(
    clin_forest$Label
  )
)

p_clin_forest <- ggplot(
  clin_forest,
  aes(
    x = HR,
    y = Label
  )
) +

  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +

  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0.16
  ) +

  geom_point(
    size = 3
  ) +

  scale_x_log10() +

  labs(
    title = "Independent prognostic value of the 8-gene risk score",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "Figure6A_Clinical_Multivariable_Cox_Forest.pdf",
  p_clin_forest,
  width = 8.5,
  height = 6
)

model_auc <- read.delim(
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalModelComparison.tsv",
  check.names = FALSE
)

auc_compare <- dplyr::bind_rows(

  data.frame(
    Model = model_auc$Model,
    Time = "1-year",
    AUC = model_auc$AUC_1y
  ),

  data.frame(
    Model = model_auc$Model,
    Time = "3-year",
    AUC = model_auc$AUC_3y
  ),

  data.frame(
    Model = model_auc$Model,
    Time = "5-year",
    AUC = model_auc$AUC_5y
  )
)

auc_compare$Time <- factor(
  auc_compare$Time,
  levels = c(
    "1-year",
    "3-year",
    "5-year"
  )
)

p_model_auc <- ggplot(
  auc_compare,
  aes(
    x = Time,
    y = AUC,
    group = Model,
    shape = Model,
    linetype = Model
  )
) +

  geom_hline(
    yintercept = 0.5,
    linetype = "dotted"
  ) +

  geom_line(
    linewidth = 0.9
  ) +

  geom_point(
    size = 3
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        AUC
      )
    ),
    vjust = -0.8,
    size = 3.3,
    show.legend = FALSE
  ) +

  coord_cartesian(
    ylim = c(
      0.58,
      0.80
    )
  ) +

  labs(
    title = "Incremental prognostic value across time",
    x = NULL,
    y = "Time-dependent AUC",
    shape = NULL,
    linetype = NULL
  ) +

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    legend.position = "bottom"
  )

ggsave(
  "Figure6C_Clinical_Molecular_Combined_AUC.pdf",
  p_model_auc,
  width = 9,
  height = 6
)

new_figs <- c(
  "Figure1D_OSCC_Hallmark_GSEA.pdf",
  "Figure2D_WGCNA_Module_OS_Cox.pdf",
  "Figure3A_Prognostic_FeatureSelection_Flow.pdf",
  "Figure3C_NestedCV_Final8_Selection.pdf",
  "Figure4A_TCGA_RiskScore_Distribution.pdf",
  "Figure4B_TCGA_SurvivalStatus.pdf",
  "Figure4E_TCGA_InternalValidation_Cindex.pdf",
  "Figure5A_GSE41613_RiskScore_Distribution.pdf",
  "Figure5B_GSE41613_SurvivalStatus.pdf",
  "Figure6A_Clinical_Multivariable_Cox_Forest.pdf",
  "Figure6C_Clinical_Molecular_Combined_AUC.pdf"
)

data.frame(
  file = new_figs,
  exists = file.exists(new_figs)
)
