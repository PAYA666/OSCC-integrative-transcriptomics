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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

deg <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_DEG_Annotated.rds"
  )
)

module_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv"
)

stopifnot(
  file.exists(module_file)
)

wgcna_gene <- read.delim(
  module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "WGCNA genes:",
  nrow(wgcna_gene),
  "\n"
)

if (
  "gene_id" %in% colnames(wgcna_gene) &&
  "gene_id" %in% colnames(deg)
) {

  wgcna_gene$gene_id_clean <- sub(
    "\\..*$",
    "",
    wgcna_gene$gene_id
  )

  deg$gene_id_clean <- sub(
    "\\..*$",
    "",
    deg$gene_id
  )

  idx_deg <- match(
    wgcna_gene$gene_id_clean,
    deg$gene_id_clean
  )

} else {

  cat(
    "gene_id unavailable; matching by gene_name.\n"
  )

  idx_deg <- match(
    wgcna_gene$gene_name,
    deg$gene_name
  )
}

wgcna_gene$padj <-
  deg$padj[idx_deg]

wgcna_gene$log2FoldChange <-
  deg$log2FoldChange[idx_deg]

wgcna_gene$Smoking_DEG <-
  !is.na(wgcna_gene$padj) &
  wgcna_gene$padj < 0.05

wgcna_gene$Strong_Smoking_DEG <-
  !is.na(wgcna_gene$padj) &
  wgcna_gene$padj < 0.05 &
  abs(wgcna_gene$log2FoldChange) >= 1

network_bg <- wgcna_gene[
  !is.na(wgcna_gene$padj),
]

cat(
  "Network genes with DEG statistics:",
  nrow(network_bg),
  "\n"
)

cat(
  "Smoking-DEGs in WGCNA:",
  sum(network_bg$Smoking_DEG),
  "\n"
)

cat(
  "Strong Smoking-DEGs in WGCNA:",
  sum(network_bg$Strong_Smoking_DEG),
  "\n"
)

modules <- unique(
  network_bg$module_color
)

enrich_list <- lapply(
  modules,
  function(mod) {

    in_module <-
      network_bg$module_color == mod

    is_deg <-
      network_bg$Smoking_DEG

    a <- sum(
      in_module &
      is_deg
    )

    b <- sum(
      in_module &
      !is_deg
    )

    c <- sum(
      !in_module &
      is_deg
    )

    d <- sum(
      !in_module &
      !is_deg
    )

    fisher_mat <- matrix(
      c(
        a, b,
        c, d
      ),
      nrow = 2,
      byrow = TRUE
    )

    ft <- fisher.test(
      fisher_mat,
      alternative = "greater"
    )

    expected_deg <-
      sum(in_module) *
      mean(is_deg)

    data.frame(

      Module = mod,

      Module_Size =
        sum(in_module),

      Smoking_DEG =
        a,

      DEG_Fraction =
        a / sum(in_module),

      Expected_DEG =
        expected_deg,

      Fold_Enrichment =
        ifelse(
          expected_deg > 0,
          a / expected_deg,
          NA_real_
        ),

      Odds_Ratio =
        unname(
          ft$estimate
        ),

      P =
        ft$p.value,

      stringsAsFactors = FALSE
    )
  }
)

module_enrich <- do.call(
  rbind,
  enrich_list
)

module_enrich$FDR <- p.adjust(
  module_enrich$P,
  method = "BH"
)

module_enrich$Significant <-
  module_enrich$FDR < 0.05

module_enrich <- module_enrich[
  order(
    module_enrich$FDR,
    -module_enrich$Odds_Ratio
  ),
]

write.table(
  module_enrich,
  file.path(
    fig_dir,
    "Figure4A_WGCNA_Module_DEG_Enrichment_SourceData_ALL.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_dat <- module_enrich[
  tolower(module_enrich$Module) != "grey",
]

plot_dat$log2_OR <- log2(
  plot_dat$Odds_Ratio
)

plot_dat$Module_Label <- paste0(
  tools::toTitleCase(
    plot_dat$Module
  ),
  "  (n=",
  plot_dat$Module_Size,
  ")"
)

plot_dat <- plot_dat[
  order(
    plot_dat$Odds_Ratio
  ),
]

plot_dat$Module_Label <- factor(
  plot_dat$Module_Label,
  levels = plot_dat$Module_Label
)

plot_dat$Plot_Group <- "Other module"

plot_dat$Plot_Group[
  plot_dat$Module == "turquoise" &
  plot_dat$FDR < 0.05
] <- "Turquoise"

plot_dat$Plot_Group[
  plot_dat$Module == "black" &
  plot_dat$FDR < 0.05
] <- "Black"

plot_dat$Plot_Group[
  plot_dat$Module == "red" &
  plot_dat$FDR < 0.05
] <- "Red"

plot_colors <- c(

  "Turquoise" =
    "#4C8D8A",

  "Black" =
    "#3F3F3F",

  "Red" =
    "#A94A44",

  "Other module" =
    "#C4C4C4"
)

format_q <- function(x) {

  if (x < 0.001) {

    paste0(
      "q=",
      formatC(
        x,
        format = "e",
        digits = 1
      )
    )

  } else {

    paste0(
      "q=",
      formatC(
        x,
        format = "f",
        digits = 3
      )
    )
  }
}

plot_dat$FDR_Label <- ifelse(
  plot_dat$FDR < 0.05,
  vapply(
    plot_dat$FDR,
    format_q,
    character(1)
  ),
  ""
)

p4a <- ggplot(
  plot_dat,
  aes(
    x = log2_OR,
    y = Module_Label
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#777777"
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = log2_OR,
      y = Module_Label,
      yend = Module_Label,
      color = Plot_Group
    ),
    linewidth = 1.15,
    alpha = 0.72
  ) +

  geom_point(
    aes(
      size = Smoking_DEG,
      fill = Plot_Group
    ),
    shape = 21,
    color = "black",
    stroke = 0.45,
    alpha = 0.96
  ) +

  geom_text(
    aes(
      label = sprintf(
        "OR %.2f",
        Odds_Ratio
      )
    ),
    hjust = ifelse(
      plot_dat$log2_OR >= 0,
      -0.55,
      1.55
    ),
    size = 3.15,
    color = "black"
  ) +

  geom_text(
    data = plot_dat[
      plot_dat$FDR < 0.05,
    ],
    aes(
      x = max(
        plot_dat$log2_OR,
        na.rm = TRUE
      ) + 0.85,
      y = Module_Label,
      label = FDR_Label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.05,
    color = "black"
  ) +

  scale_color_manual(
    values = plot_colors,
    guide = "none"
  ) +

  scale_fill_manual(
    values = plot_colors,
    guide = "none"
  ) +

  scale_size_continuous(
    range = c(
      3.4,
      7.8
    ),
    name = "Smoking-DEGs"
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.34
      )
    )
  ) +

  labs(
    x = expression(
      log[2] *
      " odds ratio for Smoking-DEG enrichment"
    ),
    y = NULL
  ) +

  theme_classic(
    base_size = 11.5,
    base_family = "sans"
  ) +

  theme(

    axis.title.x = element_text(
      size = 11.5,
      margin = margin(
        t = 8
      )
    ),

    axis.text.x = element_text(
      size = 9.5,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 10.2,
      color = "black"
    ),

    axis.line.y = element_blank(),

    axis.ticks.y = element_blank(),

    axis.line.x = element_line(
      linewidth = 0.55,
      color = "black"
    ),

    axis.ticks.x = element_line(
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
      t = 10,
      r = 55,
      b = 10,
      l = 8
    )
  ) +

  coord_cartesian(
    clip = "off"
  )

print(
  p4a
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure4A_WGCNA_SmokingDEG_Enrichment_PUBLICATION.pdf"
  ),
  plot = p4a,
  width = 8.0,
  height = 5.3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure4A_WGCNA_SmokingDEG_Enrichment_PUBLICATION.tiff"
  ),
  plot = p4a,
  width = 8.0,
  height = 5.3,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure4A_WGCNA_SmokingDEG_Enrichment_PUBLICATION.png"
  ),
  plot = p4a,
  width = 8.0,
  height = 5.3,
  units = "in",
  dpi = 200
)

cat("\n========================================\n")
cat("Figure 4A finished\n")
cat("========================================\n")

print(
  module_enrich[
    ,
    c(
      "Module",
      "Module_Size",
      "Smoking_DEG",
      "Odds_Ratio",
      "Fold_Enrichment",
      "P",
      "FDR"
    )
  ]
)
