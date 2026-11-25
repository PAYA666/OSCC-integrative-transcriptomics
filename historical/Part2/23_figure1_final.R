project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root
fig_dir  <- file.path(proj_dir, "Figures_Final")

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

pkgs <- c(
  "ggplot2",
  "ggrepel",
  "patchwork"
)

for (p in pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }

}

library(ggplot2)
library(ggrepel)
library(patchwork)

exposure <- readRDS(
  file.path(
    proj_dir,
    "Part2_TCGA_OSCC_Nature8_Exposure_341_CLEAN.rds"
  )
)

deg <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_Annotated.rds"
  )
)

n_total <- nrow(exposure)

n_never <- sum(
  exposure$smoking_bin == "Never",
  na.rm = TRUE
)

n_ever <- sum(
  exposure$smoking_bin == "Ever",
  na.rm = TRUE
)

n_missing <- sum(
  is.na(exposure$smoking_bin)
)

n_former <- sum(
  exposure$smoking3 == "Former",
  na.rm = TRUE
)

n_current <- sum(
  exposure$smoking3 == "Current",
  na.rm = TRUE
)

n_sig <- sum(
  deg$padj < 0.05,
  na.rm = TRUE
)

n_up <- sum(
  deg$padj < 0.05 &
    deg$log2FoldChange >= 1,
  na.rm = TRUE
)

n_down <- sum(
  deg$padj < 0.05 &
    deg$log2FoldChange <= -1,
  na.rm = TRUE
)

theme_pub <- theme_classic(
  base_size = 11
) +
  theme(
    text = element_text(
      family = "sans"
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9
    ),
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 9
    ),
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    legend.text = element_text(
      size = 8
    )
  )

flow_df <- data.frame(

  x = c(
    1,
    1,
    1,
    1,
    1
  ),

  y = c(
    5,
    4,
    3,
    2,
    1
  ),

  label = c(

    paste0(
      "TCGA Nature8 OSCC\n",
      n_total,
      " unique patients"
    ),

    paste0(
      "Smoking information available\n",
      n_never + n_ever,
      " patients"
    ),

    paste0(
      "Never: ",
      n_never,
      "     Ever: ",
      n_ever,
      "\n(Current ",
      n_current,
      " + Former ",
      n_former,
      ")"
    ),

    paste0(
      "Smoking-associated transcriptome\n",
      n_sig,
      " DEGs (FDR < 0.05)"
    ),

    paste0(
      "Pathway + co-expression + PPI\n",
      "integrative network analysis"
    )

  ),

  stringsAsFactors = FALSE
)

pA <- ggplot(
  flow_df,
  aes(
    x = x,
    y = y
  )
) +

  geom_segment(
    data = data.frame(
      x = 1,
      xend = 1,
      y = c(
        4.72,
        3.72,
        2.72,
        1.72
      ),
      yend = c(
        4.28,
        3.28,
        2.28,
        1.28
      )
    ),
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    arrow = arrow(
      length = unit(
        0.12,
        "inches"
      )
    )
  ) +

  geom_label(
    aes(
      label = label
    ),
    size = 3.25,
    label.size = 0.3,
    label.padding = unit(
      0.22,
      "lines"
    ),
    lineheight = 1.0,
    fill = "white"
  ) +

  coord_cartesian(
    xlim = c(
      0.25,
      1.75
    ),
    ylim = c(
      0.6,
      5.4
    )
  ) +

  labs(
    title = "A  Study design"
  ) +

  theme_void(
    base_size = 11
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    )
  )

smoking_df <- data.frame(

  Group = factor(
    c(
      "Never",
      "Ever"
    ),
    levels = c(
      "Never",
      "Ever"
    )
  ),

  N = c(
    n_never,
    n_ever
  )
)

pB <- ggplot(
  smoking_df,
  aes(
    x = Group,
    y = N,
    fill = Group
  )
) +

  geom_col(
    width = 0.62,
    linewidth = 0.4,
    color = "black"
  ) +

  geom_text(
    aes(
      label = N
    ),
    vjust = -0.5,
    size = 4,
    fontface = "bold"
  ) +

  scale_fill_manual(
    values = c(
      "Never" = "#4C78A8",
      "Ever" = "#E45756"
    )
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +

  labs(
    title = "B  Smoking exposure",
    subtitle = paste0(
      "Complete smoking information: ",
      n_never + n_ever,
      "/",
      n_total,
      " patients"
    ),
    x = NULL,
    y = "Number of patients"
  ) +

  theme_pub +

  theme(
    legend.position = "none"
  )

volcano <- deg

volcano$Volcano_group <- "Not significant"

volcano$Volcano_group[
  !is.na(volcano$padj) &
  volcano$padj < 0.05
] <- "FDR < 0.05"

volcano$Volcano_group[
  !is.na(volcano$padj) &
  volcano$padj < 0.05 &
  volcano$log2FoldChange >= 1
] <- "Up"

volcano$Volcano_group[
  !is.na(volcano$padj) &
  volcano$padj < 0.05 &
  volcano$log2FoldChange <= -1
] <- "Down"

volcano$Volcano_group <- factor(
  volcano$Volcano_group,
  levels = c(
    "Not significant",
    "FDR < 0.05",
    "Down",
    "Up"
  )
)

min_nonzero <- min(
  volcano$padj[
    volcano$padj > 0
  ],
  na.rm = TRUE
)

volcano$padj_plot <- volcano$padj

volcano$padj_plot[
  volcano$padj_plot == 0
] <- min_nonzero / 10

volcano$minus_log10_FDR <-
  -log10(
    volcano$padj_plot
  )

label_genes <- c(
  "CYP1A1",
  "CES1",
  "NTS",
  "SFRP4",
  "SFRP2",
  "PDGFRA",
  "TK1",
  "MYBL2"
)

label_df <- volcano[
  volcano$gene_name %in% label_genes,
]

pC <- ggplot(
  volcano,
  aes(
    x = log2FoldChange,
    y = minus_log10_FDR,
    color = Volcano_group
  )
) +

  geom_point(
    alpha = 0.62,
    size = 1.05
  ) +

  geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +

  geom_point(
    data = label_df,
    size = 2.1,
    color = "black"
  ) +

  geom_text_repel(
    data = label_df,
    aes(
      label = gene_name
    ),
    size = 3.0,
    color = "black",
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +

  scale_color_manual(
    values = c(
      "Not significant" = "#C7C7C7",
      "FDR < 0.05" = "#F2B447",
      "Down" = "#4C78A8",
      "Up" = "#E45756"
    )
  ) +

  labs(
    title = "C  Smoking-associated differential expression",
    subtitle = paste0(
      n_sig,
      " genes at FDR < 0.05; ",
      n_up,
      " up and ",
      n_down,
      " down at |log2FC| ≥ 1"
    ),
    x = "log2 fold change (Ever vs Never)",
    y = expression(
      -log[10](
        "FDR"
      )
    ),
    color = NULL
  ) +

  theme_pub +

  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )

fig1 <- pA |
  pB |
  pC

fig1 <- fig1 +
  plot_layout(
    widths = c(
      0.95,
      0.75,
      1.55
    )
  )

ggsave(
  filename = file.path(
    fig_dir,
    "Figure1_StudyDesign_Smoking_DEG.pdf"
  ),
  plot = fig1,
  width = 15,
  height = 5.8,
  units = "in"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure1_StudyDesign_Smoking_DEG.tiff"
  ),
  plot = fig1,
  width = 15,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  file.path(
    fig_dir,
    "Figure1A_StudyDesign.pdf"
  ),
  pA,
  width = 4.5,
  height = 5.5
)

ggsave(
  file.path(
    fig_dir,
    "Figure1B_SmokingCohort.pdf"
  ),
  pB,
  width = 4,
  height = 4.5
)

ggsave(
  file.path(
    fig_dir,
    "Figure1C_Volcano.pdf"
  ),
  pC,
  width = 7,
  height = 5.5
)

print(fig1)

cat("\n====================================\n")
cat("FINAL FIGURE 1 finished\n")
cat("====================================\n")

cat(
  "Saved to:\n",
  fig_dir,
  "\n"
)
