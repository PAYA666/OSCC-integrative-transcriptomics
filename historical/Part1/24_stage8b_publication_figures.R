rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(ggplot2)
library(dplyr)
library(survival)
library(timeROC)
library(grid)

expected_figures <- c(

  "Figure_Nature8_OSCC_WGCNA_SoftThreshold_Comparison.pdf",
  "Figure_Nature8_OSCC_WGCNA_ModuleDendrogram.pdf",
  "Figure_Nature8_OSCC_WGCNA_ModuleTraitHeatmap_FINAL.pdf",
  "Figure_Nature8_OSCC_RepeatedLASSO_Top20_PUBLICATION.pdf",

  "Figure_Final8Gene_Forest_PUBLICATION_FINAL.pdf",

  "Figure_Nature8_OSCC_Final8Gene_KM.pdf",
  "Figure_Nature8_OSCC_Final8Gene_timeROC.pdf",

  "Figure_GSE41613_Final8Gene_External_KM.pdf",
  "Figure_GSE41613_Final8Gene_External_timeROC.pdf",

  "Figure_Final8Gene_AUC_Comparison_PUBLICATION_FINAL.pdf",
  "Figure_Final8Gene_Clinical_Cindex_PUBLICATION_FINAL.pdf",

  "Figure_Nature8_OSCC_Final8Gene_Nomogram.pdf",
  "Figure_Nature8_OSCC_Final8Gene_Calibration.pdf"
)

inventory <- data.frame(
  file = expected_figures,
  exists = file.exists(expected_figures)
)

print(inventory)

workflow <- data.frame(

  step = 1:8,

  label = c(
    "TCGA-HNSC",
    "Strict anatomical filtering\nfor OSCC",
    "341 OSCC tumors\n+ 32 normal tissues",
    "DEG analysis\n+ functional enrichment",
    "WGCNA\nco-expression network",
    "Nested survival modeling\nand stability selection",
    "Frozen 8-gene\nprognostic signature",
    "GSE41613\nexternal validation"
  )
)

workflow$x <- seq_len(nrow(workflow))
workflow$y <- 1

p_workflow <- ggplot() +

  geom_segment(
    data = workflow[-nrow(workflow), ],
    aes(
      x = x + 0.32,
      xend = x + 0.68,
      y = 1,
      yend = 1
    ),
    arrow = arrow(
      length = unit(0.16, "cm")
    ),
    linewidth = 0.7
  ) +

  geom_label(
    data = workflow,
    aes(
      x = x,
      y = y,
      label = label
    ),
    size = 3.2,
    label.size = 0.4,
    label.padding = unit(0.32, "lines")
  ) +

  coord_cartesian(
    xlim = c(0.5, 8.5),
    ylim = c(0.75, 1.25),
    clip = "off"
  ) +

  labs(
    title = "Study workflow"
  ) +

  theme_void() +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    plot.margin = margin(
      20, 30, 20, 30
    )
  )

ggsave(
  "Figure1A_OSCC_StudyWorkflow.pdf",
  p_workflow,
  width = 15,
  height = 3.8
)

ggsave(
  "Figure1A_OSCC_StudyWorkflow.png",
  p_workflow,
  width = 15,
  height = 3.8,
  dpi = 600
)

deg_all <- readRDS(
  "TCGA_OSCC_Nature8_DESeq2_results.rds"
)

deg_all$PlotGroup <- "Not significant"

deg_all$PlotGroup[
  !is.na(deg_all$padj) &
  deg_all$padj < 0.05 &
  deg_all$log2FoldChange >= 1
] <- "Up"

deg_all$PlotGroup[
  !is.na(deg_all$padj) &
  deg_all$padj < 0.05 &
  deg_all$log2FoldChange <= -1
] <- "Down"

deg_all$minus_log10_padj <-
  -log10(
    pmax(
      deg_all$padj,
      1e-300,
      na.rm = TRUE
    )
  )

label_deg <- deg_all %>%
  dplyr::filter(
    PlotGroup != "Not significant",
    !is.na(gene_name)
  ) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 12)

p_volcano <- ggplot(
  deg_all,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj
  )
) +

  geom_point(
    aes(shape = PlotGroup),
    alpha = 0.45,
    size = 1
  ) +

  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +

  geom_text(
    data = label_deg,
    aes(
      label = gene_name
    ),
    size = 3,
    check_overlap = TRUE,
    vjust = -0.5
  ) +

  labs(
    title = "Differentially expressed genes in OSCC",
    subtitle = "Tumor versus normal",
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~adjusted~P),
    shape = NULL
  ) +

  theme_classic(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "top"
  )

ggsave(
  "Figure1B_OSCC_DEG_Volcano.pdf",
  p_volcano,
  width = 8,
  height = 6.5
)

ggsave(
  "Figure1B_OSCC_DEG_Volcano.png",
  p_volcano,
  width = 8,
  height = 6.5,
  dpi = 600
)

library(clusterProfiler)
library(msigdbr)

rank_df <- deg_all[
  !is.na(deg_all$log2FoldChange) &
  !is.na(deg_all$gene_name) &
  deg_all$gene_name != "",
  c(
    "gene_name",
    "log2FoldChange"
  )
]

rank_df <- rank_df[
  order(
    abs(rank_df$log2FoldChange),
    decreasing = TRUE
  ),
]

rank_df <- rank_df[
  !duplicated(rank_df$gene_name),
]

gene_list <- rank_df$log2FoldChange
names(gene_list) <- rank_df$gene_name
gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

hallmark <- msigdbr::msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_term2gene <- hallmark[
  ,
  c(
    "gs_name",
    "gene_symbol"
  )
]

set.seed(20260813)

gsea_h <- clusterProfiler::GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_term2gene,
  pvalueCutoff = 1,
  verbose = FALSE
)

gsea_df <- as.data.frame(gsea_h)

gsea_sig <- gsea_df %>%
  dplyr::filter(
    p.adjust < 0.05
  ) %>%
  dplyr::arrange(p.adjust)

gsea_top <- gsea_sig %>%
  dplyr::slice_head(
    n = min(10, nrow(gsea_sig))
  )

gsea_top$Description <- gsub(
  "^HALLMARK_",
  "",
  gsea_top$Description
)

gsea_top$Description <- gsub(
  "_",
  " ",
  gsea_top$Description
)

gsea_top$Description <- factor(
  gsea_top$Description,
  levels = rev(gsea_top$Description)
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
    linewidth = 0.6
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = NES,
      y = Description,
      yend = Description
    ),
    linewidth = 1
  ) +

  geom_point(
    size = 3
  ) +

  labs(
    title = "Hallmark pathway enrichment in OSCC",
    x = "Normalized enrichment score (NES)",
    y = NULL
  ) +

  theme_classic(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 10
    )
  )

ggsave(
  "Figure1D_OSCC_Hallmark_GSEA.pdf",
  p_gsea,
  width = 9,
  height = 6.5
)

ggsave(
  "Figure1D_OSCC_Hallmark_GSEA.png",
  p_gsea,
  width = 9,
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

module_cox$module <- factor(
  module_cox$module,
  levels = rev(
    module_cox$module
  )
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
    linetype = "dashed"
  ) +

  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0.18
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

  theme_classic(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  "Figure2D_WGCNA_Module_OS_Cox.pdf",
  p_module_cox,
  width = 7.5,
  height = 5.8
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

n_deg <- nrow(deg_df)

n_non_grey <- sum(
  gene_module$module_color != "grey"
)

candidate_non_grey_deg <- gene_module %>%
  dplyr::filter(
    module_color != "grey",
    gene_id %in% deg_df$gene_id
  ) %>%
  dplyr::distinct(gene_id)

n_candidate <- nrow(
  candidate_non_grey_deg
)

n_stable8 <- sum(
  nested_selection$selection_frequency >= 0.80
)

funnel <- data.frame(
  x = 1:4,
  y = 1,
  label = c(
    paste0(
      "DEGs\nn = ",
      n_deg
    ),
    paste0(
      "Non-grey WGCNA genes\nn = ",
      n_non_grey
    ),
    paste0(
      "DEG ∩ non-grey WGCNA\nn = ",
      n_candidate
    ),
    paste0(
      "Stable final signature\nn = ",
      n_stable8
    )
  )
)

p_funnel <- ggplot() +

  geom_segment(
    data = funnel[-nrow(funnel), ],
    aes(
      x = x + 0.30,
      xend = x + 0.70,
      y = 1,
      yend = 1
    ),
    arrow = arrow(
      length = unit(0.18, "cm")
    ),
    linewidth = 0.8
  ) +

  geom_label(
    data = funnel,
    aes(
      x = x,
      y = y,
      label = label
    ),
    size = 4,
    label.padding = unit(
      0.4,
      "lines"
    )
  ) +

  coord_cartesian(
    xlim = c(
      0.5,
      4.5
    ),
    ylim = c(
      0.75,
      1.25
    ),
    clip = "off"
  ) +

  labs(
    title = "Feature-selection strategy"
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
  p_funnel,
  width = 11,
  height = 3.5
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
    linetype = "dashed"
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

tcga <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_RiskScore.rds"
)

tcga_plot <- tcga %>%
  dplyr::arrange(
    risk_score
  )

tcga_plot$rank <- seq_len(
  nrow(tcga_plot)
)

p_tcga_risk <- ggplot(
  tcga_plot,
  aes(
    x = rank,
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
    x = rank,
    y = OS_time
  )
) +

  geom_point(
    aes(
      shape = factor(OS_event)
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
  Evaluation = c(
    "TCGA apparent",
    "TCGA nested-CV (OOF)"
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

geo_plot$rank <- seq_len(
  nrow(geo_plot)
)

p_geo_risk <- ggplot(
  geo_plot,
  aes(
    x = rank,
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
    x = rank,
    y = OS_time_days
  )
) +

  geom_point(
    aes(
      shape = factor(OS_event)
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

clin_forest$Variable <- factor(
  clin_forest$Variable,
  levels = rev(
    clin_forest$Variable
  )
)

p_clin_forest <- ggplot(
  clin_forest,
  aes(
    x = HR,
    y = Variable
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
  width = 8,
  height = 6
)

model_auc <- read.delim(
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalModelComparison.tsv",
  check.names = FALSE
)

auc_compare <- bind_rows(

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

  coord_cartesian(
    ylim = c(
      0.58,
      0.80
    )
  ) +

  labs(
    title = "Added prognostic value of the molecular signature",
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

cat("\n========================================\n")
cat("NEW MAIN FIGURES CREATED\n")
cat("========================================\n")

new_figs <- c(
  "Figure1A_OSCC_StudyWorkflow.pdf",
  "Figure1B_OSCC_DEG_Volcano.pdf",
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

print(
  data.frame(
    file = new_figs,
    exists = file.exists(new_figs)
  )
)
