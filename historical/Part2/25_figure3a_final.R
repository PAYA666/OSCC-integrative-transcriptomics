project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root
fig_dir  <- file.path(proj_dir, "Figures_Final")

dir.create(
  fig_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

extract_gsea_df <- function(x) {

  if (inherits(x, "gseaResult")) {
    x <- as.data.frame(x)
  }

  if (is.data.frame(x)) {
    return(x)
  }

  return(NULL)
}

gsea_df <- NULL
gsea_source <- NULL

possible_objects <- c(
  "gsea_res",
  "gsea_hallmark",
  "gsea_result",
  "gsea_results",
  "GSEA_res",
  "GSEA_result"
)

for (obj_name in possible_objects) {

  if (exists(obj_name, inherits = TRUE)) {

    tmp <- try(
      extract_gsea_df(
        get(obj_name, inherits = TRUE)
      ),
      silent = TRUE
    )

    if (
      !inherits(tmp, "try-error") &&
      !is.null(tmp) &&
      "NES" %in% colnames(tmp)
    ) {

      gsea_df <- tmp
      gsea_source <- paste0(
        "workspace object: ",
        obj_name
      )

      break
    }
  }
}

if (is.null(gsea_df)) {

  gsea_files <- list.files(
    proj_dir,
    pattern = "GSEA|Hallmark",
    full.names = TRUE,
    ignore.case = TRUE
  )

  priority_files <- gsea_files[
    grepl(
      "GSEA",
      basename(gsea_files),
      ignore.case = TRUE
    ) &
    grepl(
      "Hallmark",
      basename(gsea_files),
      ignore.case = TRUE
    )
  ]

  search_files <- unique(
    c(
      priority_files,
      gsea_files
    )
  )

  for (f in search_files) {

    ext <- tolower(
      tools::file_ext(f)
    )

    tmp <- try({

      if (ext == "rds") {

        extract_gsea_df(
          readRDS(f)
        )

      } else if (ext %in% c("tsv", "txt")) {

        read.delim(
          f,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

      } else if (ext == "csv") {

        read.csv(
          f,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

      } else {

        NULL
      }

    }, silent = TRUE)

    if (
      !inherits(tmp, "try-error") &&
      !is.null(tmp) &&
      is.data.frame(tmp) &&
      "NES" %in% colnames(tmp)
    ) {

      gsea_df <- tmp
      gsea_source <- f
      break
    }
  }
}

if (is.null(gsea_df)) {

  cat("\nPossible GSEA/Hallmark files:\n")

  print(
    list.files(
      proj_dir,
      pattern = "GSEA|Hallmark",
      ignore.case = TRUE
    )
  )

  stop(
    paste0(
      "\nCould not automatically locate the Hallmark GSEA table.\n",
      "The result must contain an NES column."
    )
  )
}

cat(
  "\nGSEA source:\n",
  gsea_source,
  "\n"
)

pathway_candidates <- c(
  "Description",
  "ID",
  "Pathway",
  "pathway",
  "Term",
  "term"
)

pathway_col <- pathway_candidates[
  pathway_candidates %in%
    colnames(gsea_df)
][1]

fdr_candidates <- c(
  "p.adjust",
  "padj",
  "FDR",
  "qvalue",
  "q.value"
)

fdr_col <- fdr_candidates[
  fdr_candidates %in%
    colnames(gsea_df)
][1]

if (
  is.na(pathway_col) ||
  is.na(fdr_col)
) {

  cat("\nGSEA columns:\n")
  print(colnames(gsea_df))

  stop(
    "Could not identify pathway-name or FDR column."
  )
}

plot_gsea <- data.frame(

  Pathway = as.character(
    gsea_df[[pathway_col]]
  ),

  NES = as.numeric(
    gsea_df$NES
  ),

  FDR = as.numeric(
    gsea_df[[fdr_col]]
  ),

  stringsAsFactors = FALSE
)

plot_gsea <- plot_gsea[
  !is.na(plot_gsea$NES) &
  !is.na(plot_gsea$FDR),
]

plot_gsea <- plot_gsea[
  plot_gsea$FDR < 0.05,
]

if (nrow(plot_gsea) == 0) {
  stop("No Hallmark pathway has FDR < 0.05.")
}

plot_gsea$Pathway_clean <- gsub(
  "^HALLMARK_",
  "",
  plot_gsea$Pathway
)

plot_gsea$Pathway_clean <- gsub(
  "_",
  " ",
  plot_gsea$Pathway_clean
)

cellcycle_terms <- c(
  "E2F TARGETS",
  "G2M CHECKPOINT",
  "MYC TARGETS V1",
  "MYC TARGETS V2",
  "MITOTIC SPINDLE"
)

metabolic_terms <- c(
  "OXIDATIVE PHOSPHORYLATION",
  "XENOBIOTIC METABOLISM",
  "FATTY ACID METABOLISM",
  "GLYCOLYSIS",
  "ADIPOGENESIS"
)

plot_gsea$Program <- "Signaling / other"

plot_gsea$Program[
  toupper(plot_gsea$Pathway_clean) %in%
    cellcycle_terms
] <- "Cell cycle / proliferation"

plot_gsea$Program[
  toupper(plot_gsea$Pathway_clean) %in%
    metabolic_terms
] <- "Metabolism / exposure"

plot_gsea$Program <- factor(
  plot_gsea$Program,
  levels = c(
    "Cell cycle / proliferation",
    "Metabolism / exposure",
    "Signaling / other"
  )
)

key_terms <- c(
  "E2F TARGETS",
  "G2M CHECKPOINT",
  "MYC TARGETS V1",
  "OXIDATIVE PHOSPHORYLATION",
  "XENOBIOTIC METABOLISM",
  "MTORC1 SIGNALING"
)

plot_gsea$Key <- toupper(
  plot_gsea$Pathway_clean
) %in% key_terms

plot_gsea$minus_log10_FDR <- -log10(
  pmax(
    plot_gsea$FDR,
    1e-300
  )
)

plot_gsea$SizeScore <- pmin(
  plot_gsea$minus_log10_FDR,
  12
)

plot_gsea <- plot_gsea[
  order(
    plot_gsea$NES
  ),
]

plot_gsea$Pathway_clean <- factor(
  plot_gsea$Pathway_clean,
  levels = plot_gsea$Pathway_clean
)

program_colors <- c(

  "Cell cycle / proliferation" =
    "#A94A44",

  "Metabolism / exposure" =
    "#4E7A70",

  "Signaling / other" =
    "#8A8A8A"
)

p_gsea <- ggplot(
  plot_gsea,
  aes(
    x = NES,
    y = Pathway_clean
  )
) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    color = "#8A8A8A",
    linetype = "dashed"
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = NES,
      y = Pathway_clean,
      yend = Pathway_clean,
      color = Program
    ),
    linewidth = 1.15,
    alpha = 0.65
  ) +

  geom_point(
    aes(
      size = SizeScore,
      fill = Program
    ),
    shape = 21,
    color = "white",
    stroke = 0.45,
    alpha = 0.95
  ) +

  geom_point(
    data = plot_gsea[
      plot_gsea$Key,
    ],
    aes(
      x = NES,
      y = Pathway_clean,
      size = SizeScore
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 0.9
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        NES
      )
    ),
    hjust = ifelse(
      plot_gsea$NES >= 0,
      -0.55,
      1.55
    ),
    size = 3.0,
    color = "black"
  ) +

  scale_color_manual(
    values = program_colors
  ) +

  scale_fill_manual(
    values = program_colors
  ) +

  scale_size_continuous(
    range = c(
      3.2,
      7.8
    ),
    name = expression(
      -log[10] *
        "(FDR)"
    )
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(
        0.03,
        0.13
      )
    )
  ) +

  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    fill = NULL,
    color = NULL
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
      size = 9.5,
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

    legend.text = element_text(
      size = 9
    ),

    legend.title = element_text(
      size = 9.5
    ),

    plot.margin = margin(
      t = 10,
      r = 25,
      b = 10,
      l = 8
    )
  )

print(
  p_gsea
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure3A_HallmarkGSEA_PUBLICATION.pdf"
  ),
  plot = p_gsea,
  width = 8.0,
  height = 6.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure3A_HallmarkGSEA_PUBLICATION.tiff"
  ),
  plot = p_gsea,
  width = 8.0,
  height = 6.6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

if (
  requireNamespace(
    "svglite",
    quietly = TRUE
  )
) {

  ggsave(
    filename = file.path(
      fig_dir,
      "Figure3A_HallmarkGSEA_PUBLICATION.svg"
    ),
    plot = p_gsea,
    width = 8.0,
    height = 6.6,
    units = "in"
  )
}

fig3a_source <- plot_gsea[
  ,
  c(
    "Pathway",
    "Pathway_clean",
    "NES",
    "FDR",
    "Program",
    "Key"
  )
]

write.table(
  fig3a_source,
  file.path(
    fig_dir,
    "Figure3A_SourceData.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Figure 3A Hallmark GSEA finished\n")
cat("========================================\n")

cat(
  "Significant Hallmark pathways:",
  nrow(plot_gsea),
  "\n"
)

cat(
  "Ever-enriched:",
  sum(
    plot_gsea$NES > 0
  ),
  "\n"
)

cat(
  "Never-enriched:",
  sum(
    plot_gsea$NES < 0
  ),
  "\n"
)

print(
  fig3a_source[
    order(
      -fig3a_source$NES
    ),
  ]
)
