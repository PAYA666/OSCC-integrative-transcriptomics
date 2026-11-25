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

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}

library(cowplot)

gsea_main <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_Hallmark_GSEA_MAIN.rds"
  )
)

gsea_df <- as.data.frame(
  gsea_main
)

pathway_show <- c(

  "HALLMARK_XENOBIOTIC_METABOLISM",

  "HALLMARK_E2F_TARGETS",

  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"

)

short_names <- c(

  "Xenobiotic metabolism",

  "E2F targets",

  "Oxidative phosphorylation"

)

missing_pathways <- setdiff(
  pathway_show,
  gsea_df$ID
)

if (length(missing_pathways) > 0) {

  cat("Missing pathways:\n")

  print(
    missing_pathways
  )

  stop(
    "Requested pathway missing from Smoking MAIN GSEA."
  )
}

title_info <- gsea_df[
  match(
    pathway_show,
    gsea_df$ID
  ),
  c(
    "ID",
    "Description",
    "NES",
    "p.adjust"
  )
]

title_info$ShortName <-
  short_names

format_fdr <- function(x) {

  if (x < 0.001) {

    formatC(
      x,
      format = "e",
      digits = 1
    )

  } else {

    formatC(
      x,
      format = "f",
      digits = 3
    )

  }
}

make_gsea_panel <- function(
  gs_id,
  short_name,
  nes,
  fdr,
  curve_color
) {

  title_text <- paste0(
    short_name,
    "   |   NES = ",
    sprintf(
      "%.2f",
      nes
    ),
    ", FDR = ",
    format_fdr(
      fdr
    )
  )

  p <- enrichplot::gseaplot2(

    gsea_main,

    geneSetID = gs_id,

    title = title_text,

    color = curve_color,

    base_size = 11,

    rel_heights = c(
      1.35,
      0.35,
      0.80
    ),

    subplots = 1:3,

    pvalue_table = FALSE

  )

  return(
    p
  )
}

col_xeno  <- "#4E7A70"
col_e2f   <- "#A94A44"
col_oxphos <- "#596F8A"

p_xeno <- make_gsea_panel(

  gs_id =
    title_info$ID[1],

  short_name =
    title_info$ShortName[1],

  nes =
    title_info$NES[1],

  fdr =
    title_info$p.adjust[1],

  curve_color =
    col_xeno
)

p_e2f <- make_gsea_panel(

  gs_id =
    title_info$ID[2],

  short_name =
    title_info$ShortName[2],

  nes =
    title_info$NES[2],

  fdr =
    title_info$p.adjust[2],

  curve_color =
    col_e2f
)

p_oxphos <- make_gsea_panel(

  gs_id =
    title_info$ID[3],

  short_name =
    title_info$ShortName[3],

  nes =
    title_info$NES[3],

  fdr =
    title_info$p.adjust[3],

  curve_color =
    col_oxphos
)

print(
  p_xeno
)

print(
  p_e2f
)

print(
  p_oxphos
)

fig3b <- cowplot::plot_grid(

  p_xeno,

  p_e2f,

  p_oxphos,

  ncol = 1,

  labels = c(
    "A",
    "B",
    "C"
  ),

  label_size = 13,

  label_fontface = "bold",

  align = "v",

  axis = "lr",

  rel_heights = c(
    1,
    1,
    1
  )

)

print(
  fig3b
)

ggsave(

  filename = file.path(
    fig_dir,
    "Figure3B_Representative_GSEA_Curves_PUBLICATION.pdf"
  ),

  plot = fig3b,

  width = 8.2,

  height = 11.2,

  units = "in",

  device = cairo_pdf

)

ggsave(

  filename = file.path(
    fig_dir,
    "Figure3B_Representative_GSEA_Curves_PUBLICATION.tiff"
  ),

  plot = fig3b,

  width = 8.2,

  height = 11.2,

  units = "in",

  dpi = 600,

  compression = "lzw"

)

ggsave(

  file.path(
    fig_dir,
    "Figure3B1_Xenobiotic_GSEA.pdf"
  ),

  p_xeno,

  width = 8.0,

  height = 3.7,

  units = "in",

  device = cairo_pdf

)

ggsave(

  file.path(
    fig_dir,
    "Figure3B2_E2F_GSEA.pdf"
  ),

  p_e2f,

  width = 8.0,

  height = 3.7,

  units = "in",

  device = cairo_pdf

)

ggsave(

  file.path(
    fig_dir,
    "Figure3B3_OXPHOS_GSEA.pdf"
  ),

  p_oxphos,

  width = 8.0,

  height = 3.7,

  units = "in",

  device = cairo_pdf

)

write.table(

  title_info,

  file.path(
    fig_dir,
    "Figure3B_SourceData.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

cat("\n========================================\n")
cat("Figure 3B finished\n")
cat("========================================\n")

print(
  title_info[
    ,
    c(
      "ShortName",
      "NES",
      "p.adjust"
    )
  ]
)
