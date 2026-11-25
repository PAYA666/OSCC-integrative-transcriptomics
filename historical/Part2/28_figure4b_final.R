project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root
fig_dir  <- file.path(
  proj_dir,
  "Figures_Final"
)

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

pkgs <- c(
  "ggplot2",
  "ggrepel"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(ggplot2)
library(ggrepel)

cand <- data.frame(

  gene_name = c(
    "SFRP4",
    "PDGFRA",
    "SFRP2",
    "TK1",
    "OGN",
    "MYBL2",
    "FNDC1"
  ),

  Module = c(
    "turquoise",
    "turquoise",
    "turquoise",
    "black",
    "turquoise",
    "black",
    "turquoise"
  ),

  log2FoldChange = c(
    1.2391988,
    0.5142705,
    0.6660361,
    0.3539956,
    1.1177555,
    0.3430908,
    0.7861249
  ),

  padj = c(
    1.223027e-03,
    4.552026e-02,
    3.019642e-02,
    1.516732e-02,
    4.534392e-03,
    1.790387e-02,
    1.624418e-02
  ),

  MM = c(
    0.761,
    0.826,
    0.837,
    0.725,
    0.769,
    0.716,
    0.838
  ),

  GS = c(
    0.103,
    0.138,
    0.108,
    0.167,
    0.106,
    0.170,
    0.117
  ),

  PPI_Degree = c(
    14,
    9,
    10,
    6,
    9,
    5,
    6
  ),

  MeanTopologyRank = c(
    7.67,
    9.00,
    12.00,
    13.33,
    15.67,
    18.00,
    20.00
  ),

  stringsAsFactors = FALSE
)

cand$abs_MM <- abs(
  cand$MM
)

cand$minus_log10_FDR <- -log10(
  cand$padj
)

module_colors <- c(

  "turquoise" =
    "#4C8D8A",

  "black" =
    "#444444"
)

cand$label <- cand$gene_name

p4b <- ggplot(
  cand,
  aes(
    x = abs_MM,
    y = PPI_Degree
  )
) +

  geom_vline(
    xintercept = 0.70,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#8A8A8A"
  ) +

  geom_hline(
    yintercept = c(
      5,
      10,
      15
    ),
    linewidth = 0.3,
    color = "#E5E5E5"
  ) +

  geom_point(
    aes(
      size = minus_log10_FDR,
      fill = Module
    ),
    shape = 21,
    color = "black",
    stroke = 0.65,
    alpha = 0.96
  ) +

  geom_text_repel(
    aes(
      label = label
    ),
    size = 3.7,
    fontface = "italic",
    color = "black",

    box.padding = 0.55,
    point.padding = 0.45,

    segment.color = "#777777",
    segment.size = 0.35,

    min.segment.length = 0,
    max.overlaps = Inf,

    seed = 123
  ) +

  scale_fill_manual(
    values = module_colors,
    name = "WGCNA module",
    labels = c(
      "black" = "Black",
      "turquoise" = "Turquoise"
    )
  ) +

  scale_size_continuous(
    range = c(
      4.5,
      8.5
    ),
    name = expression(
      -log[10] *
        "(DEG FDR)"
    )
  ) +

  scale_x_continuous(
    limits = c(
      0.69,
      0.86
    ),
    breaks = seq(
      0.70,
      0.85,
      by = 0.05
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.05
      )
    )
  ) +

  scale_y_continuous(
    limits = c(
      3.5,
      15.5
    ),
    breaks = c(
      5,
      10,
      15
    ),
    expand = expansion(
      mult = c(
        0.03,
        0.05
      )
    )
  ) +

  labs(
    x = expression(
      "|" * MM * "|  (module membership)"
    ),
    y = "PPI degree"
  ) +

  theme_classic(
    base_size = 11.5,
    base_family = "sans"
  ) +

  theme(

    axis.title = element_text(
      size = 11.5,
      color = "black"
    ),

    axis.text = element_text(
      size = 9.5,
      color = "black"
    ),

    axis.line = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),

    legend.position = "right",

    legend.title = element_text(
      size = 9.5
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      12,
      25,
      10,
      10
    )
  ) +

  coord_cartesian(
    clip = "off"
  )

print(
  p4b
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure4B_MultiEvidence_Candidates_PUBLICATION.pdf"
  ),
  plot = p4b,
  width = 7.0,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure4B_MultiEvidence_Candidates_PUBLICATION.tiff"
  ),
  plot = p4b,
  width = 7.0,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure4B_MultiEvidence_Candidates_PUBLICATION.png"
  ),
  plot = p4b,
  width = 7.0,
  height = 5.8,
  units = "in",
  dpi = 200
)

write.table(
  cand,
  file.path(
    fig_dir,
    "Figure4B_SourceData.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Figure 4B finished\n")
cat("========================================\n")

print(
  cand[
    order(
      -cand$PPI_Degree
    ),
    c(
      "gene_name",
      "Module",
      "log2FoldChange",
      "padj",
      "MM",
      "GS",
      "PPI_Degree",
      "MeanTopologyRank"
    )
  ]
)
