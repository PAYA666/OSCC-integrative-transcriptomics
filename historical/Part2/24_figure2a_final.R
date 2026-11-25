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
  "patchwork"
)

for (p in pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }

}

library(ggplot2)
library(patchwork)

expr_long <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_VST_LongFormat.rds"
  )
)

exposure <- readRDS(
  file.path(
    proj_dir,
    "Part2_TCGA_OSCC_Nature8_Exposure_341_CLEAN.rds"
  )
)

if (!"smoking3" %in% colnames(exposure)) {

  status_raw <- trimws(
    as.character(
      exposure$tobacco_smoking_status
    )
  )

  status_lower <- tolower(
    status_raw
  )

  exposure$smoking3 <- NA_character_

  exposure$smoking3[
    !is.na(status_lower) &
      (
        grepl("lifelong", status_lower) |
        grepl("never", status_lower)
      )
  ] <- "Never"

  exposure$smoking3[
    !is.na(status_lower) &
      grepl("current", status_lower)
  ] <- "Current"

  exposure$smoking3[
    !is.na(status_lower) &
      (
        grepl("former", status_lower) |
        grepl("reformed", status_lower)
      )
  ] <- "Former"
}

exposure$smoking3 <- factor(
  exposure$smoking3,
  levels = c(
    "Never",
    "Former",
    "Current"
  )
)

meta_sub <- exposure[
  ,
  c(
    "patient_id",
    "smoking3"
  )
]

plot_dat <- merge(
  expr_long,
  meta_sub,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

plot_dat <- plot_dat[
  !is.na(plot_dat$smoking3),
]

sens <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_ExposureDose_Sensitivity_FINAL.rds"
  )
)

pairwise_results <- sens$pairwise_results

exposure_genes <- c(
  "CYP1A1",
  "CES1"
)

prolif_genes <- c(
  "TK1",
  "MYBL2"
)

selected_genes <- c(
  exposure_genes,
  prolif_genes
)

plot_dat <- plot_dat[
  plot_dat$gene_name %in% selected_genes,
]

plot_dat$gene_name <- factor(
  plot_dat$gene_name,
  levels = selected_genes
)

group_n <- table(
  exposure$smoking3
)

x_labels <- c(

  "Never" = paste0(
    "Never\nn = ",
    group_n["Never"]
  ),

  "Former" = paste0(
    "Former\nn = ",
    group_n["Former"]
  ),

  "Current" = paste0(
    "Current\nn = ",
    group_n["Current"]
  )
)

group_colors <- c(

  "Never"   = "#B8B8B8",

  "Former"  = "#C99A4A",

  "Current" = "#B84C4C"
)

format_q <- function(x) {

  if (x < 0.001) {

    return(
      paste0(
        "q=",
        formatC(
          x,
          format = "e",
          digits = 1
        )
      )
    )

  } else {

    return(
      paste0(
        "q=",
        formatC(
          x,
          format = "f",
          digits = 3
        )
      )
    )

  }
}

make_annotations <- function(genes) {

  sig <- pairwise_results[
    pairwise_results$gene_name %in% genes &
      pairwise_results$FDR < 0.05,
  ]

  if (nrow(sig) == 0) {
    return(NULL)
  }

  group_position <- c(
    "Never" = 1,
    "Former" = 2,
    "Current" = 3
  )

  ann_list <- list()

  for (g in genes) {

    sig_g <- sig[
      sig$gene_name == g,
    ]

    if (nrow(sig_g) == 0) {
      next
    }

    dat_g <- plot_dat[
      plot_dat$gene_name == g,
    ]

    y_min <- min(
      dat_g$VST_expression,
      na.rm = TRUE
    )

    y_max <- max(
      dat_g$VST_expression,
      na.rm = TRUE
    )

    y_range <- y_max - y_min

    step <- max(
      y_range * 0.10,
      0.20
    )

    desired_order <- c(
      "Former_vs_Never",
      "Current_vs_Never",
      "Current_vs_Former"
    )

    sig_g$ord <- match(
      sig_g$Comparison,
      desired_order
    )

    sig_g <- sig_g[
      order(sig_g$ord),
    ]

    for (i in seq_len(nrow(sig_g))) {

      comp <- strsplit(
        sig_g$Comparison[i],
        "_vs_"
      )[[1]]

      group1 <- comp[1]
      group2 <- comp[2]

      x1 <- group_position[group1]
      x2 <- group_position[group2]

      ann_list[[
        length(ann_list) + 1
      ]] <- data.frame(

        gene_name = g,

        x1 = min(
          x1,
          x2
        ),

        x2 = max(
          x1,
          x2
        ),

        y =
          y_max +
          step * i,

        label =
          format_q(
            sig_g$FDR[i]
          ),

        stringsAsFactors = FALSE
      )

    }

  }

  do.call(
    rbind,
    ann_list
  )
}

make_gene_row <- function(
  genes,
  row_title
) {

  dat <- plot_dat[
    plot_dat$gene_name %in% genes,
  ]

  dat$gene_name <- factor(
    dat$gene_name,
    levels = genes
  )

  ann <- make_annotations(
    genes
  )

  p <- ggplot(
    dat,
    aes(
      x = smoking3,
      y = VST_expression
    )
  ) +

    geom_jitter(
      aes(
        color = smoking3
      ),
      width = 0.13,
      height = 0,
      size = 1.05,
      alpha = 0.32,
      show.legend = FALSE
    ) +

    geom_boxplot(
      aes(
        fill = smoking3
      ),
      width = 0.50,
      outlier.shape = NA,
      alpha = 0.72,
      linewidth = 0.55,
      color = "black",
      show.legend = FALSE
    ) +

    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 2.3,
      stroke = 0.6,
      fill = "white",
      color = "black"
    ) +

    facet_wrap(
      ~ gene_name,
      nrow = 1,
      scales = "free_y"
    ) +

    scale_fill_manual(
      values = group_colors
    ) +

    scale_color_manual(
      values = group_colors
    ) +

    scale_x_discrete(
      labels = x_labels
    ) +

    labs(
      title = row_title,
      x = NULL,
      y = "VST expression"
    ) +

    theme_classic(
      base_size = 11,
      base_family = "sans"
    ) +

    theme(

      plot.title = element_text(
        size = 11.5,
        face = "bold",
        hjust = 0,
        margin = margin(
          b = 7
        )
      ),

      strip.background = element_blank(),

      strip.text = element_text(
        size = 11.5,
        face = "bold.italic",
        color = "black",
        margin = margin(
          b = 4
        )
      ),

      axis.title.y = element_text(
        size = 10.5
      ),

      axis.text.x = element_text(
        size = 9,
        color = "black",
        lineheight = 0.9
      ),

      axis.text.y = element_text(
        size = 9,
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

      panel.spacing = unit(
        1.1,
        "lines"
      ),

      plot.margin = margin(
        8,
        12,
        8,
        8
      )
    )

  if (!is.null(ann)) {

    ann$gene_name <- factor(
      ann$gene_name,
      levels = genes
    )

    tick_height <- 0.04

    p <- p +

      geom_segment(
        data = ann,
        aes(
          x = x1,
          xend = x2,
          y = y,
          yend = y
        ),
        inherit.aes = FALSE,
        linewidth = 0.45,
        color = "black"
      ) +

      geom_segment(
        data = ann,
        aes(
          x = x1,
          xend = x1,
          y = y,
          yend = y - tick_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.45,
        color = "black"
      ) +

      geom_segment(
        data = ann,
        aes(
          x = x2,
          xend = x2,
          y = y,
          yend = y - tick_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.45,
        color = "black"
      ) +

      geom_text(
        data = ann,
        aes(
          x = (
            x1 + x2
          ) / 2,
          y = y,
          label = label
        ),
        inherit.aes = FALSE,
        vjust = -0.40,
        size = 3.0,
        color = "black"
      )

  }

  p <- p +
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.04,
          0.20
        )
      )
    )

  p
}

p_exposure <- make_gene_row(

  genes = exposure_genes,

  row_title =
    "Exposure-response genes"
)

p_prolif <- make_gene_row(

  genes = prolif_genes,

  row_title =
    "Proliferative program"
)

fig2A <- p_exposure /
  p_prolif +
  plot_layout(
    heights = c(
      1,
      1
    )
  )

print(
  fig2A
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure2A_SmokingStatus_KeyGenes_PUBLICATION.pdf"
  ),
  plot = fig2A,
  width = 7.2,
  height = 7.2,
  units = "in"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure2A_SmokingStatus_KeyGenes_PUBLICATION.tiff"
  ),
  plot = fig2A,
  width = 7.2,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure2A1_ExposureResponse_PUBLICATION.pdf"
  ),
  plot = p_exposure,
  width = 7.2,
  height = 3.5,
  units = "in"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure2A2_ProliferativeProgram_PUBLICATION.pdf"
  ),
  plot = p_prolif,
  width = 7.2,
  height = 3.5,
  units = "in"
)

cat("\n========================================\n")
cat("Figure 2A finished\n")
cat("========================================\n")

cat(
  "Saved to:",
  fig_dir,
  "\n"
)
