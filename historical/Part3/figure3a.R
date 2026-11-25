project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure3")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ultra7 <- c(
"CCN1",
"CCL2",
"ZBTB10",
"DAB2",
"JUN",
"CCL20",
"KYNU"
)

plot_df <- step4$candidate_module_mapping %>%

  filter(
    Gene %in% ultra7
  ) %>%

  distinct(
    Gene,
    module_color
  )

plot_df$Gene <- factor(
  plot_df$Gene,
  levels = rev(ultra7)
)

print(plot_df)

p <- ggplot(
  plot_df,
  aes(
    x=module_color,
    y=Gene
  )
)+

geom_point(
  size=8
)+

theme_classic(
  base_size=14
)+

labs(
  title=
  "TCGA WGCNA projection of frozen Ultra7 program",

  x=
  "WGCNA module",

  y=NULL
)

ggsave(
  file.path(
    out_dir,
    "Figure3A_Ultra7_TCGA_WGCNA_projection.png"
  ),

  p,

  width=90,
  height=65,
  units="mm",
  dpi=600,
  bg="white"
)

print(p)
