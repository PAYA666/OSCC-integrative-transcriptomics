project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure3")

ultra_network <- step5$ultra_network %>%

  filter(
    Gene %in% c(
      "CCN1",
      "CCL2",
      "ZBTB10",
      "DAB2",
      "JUN",
      "CCL20",
      "KYNU"
    )
  )

evidence_df <- ultra_network %>%

  transmute(

    Gene,

    Meta_Z,

    Cross_model_DEG =
      as.numeric(Ultra_Strict),

    Pathway_convergence =
      BothDataset_Pathway_Count,

    TCGA_WGCNA =
      as.numeric(
        In_TCGA_VST &
        In_WGCNA_top5000
      ),

    Network_support =
      as.numeric(
        Network_Interpretable
      )
  )

plot_df <- evidence_df %>%

  pivot_longer(
    cols=
      c(
        Cross_model_DEG,
        Pathway_convergence,
        TCGA_WGCNA,
        Network_support
      ),

    names_to="Evidence",

    values_to="Score"
  )

plot_df$Gene <-
factor(
  plot_df$Gene,
  levels=c(
    "CCN1",
    "CCL2",
    "ZBTB10",
    "DAB2",
    "JUN",
    "CCL20",
    "KYNU"
  )
)

p <- ggplot(
  plot_df,
  aes(
    x=Evidence,
    y=Gene
  )
)+

geom_point(
  aes(
    size=Score,
    color=Meta_Z
  )
)+

scale_color_gradient(
  low="grey70",
  high="firebrick"
)+

theme_classic(
  base_size=14
)+

theme(
  axis.text.x =
    element_text(
      angle=45,
      hjust=1
    )
)+

labs(
  title=
    "Multi-layer evidence supporting frozen Ultra7 candidates",

  x=NULL,

  y=NULL,

  size="Evidence",

  color="Meta Z"
)

ggsave(
  file.path(
    out_dir,
    "Figure3C_Ultra7_multilayer_support.png"
  ),

  p,

  width=95,
  height=70,
  units="mm",
  dpi=600,
  bg="white"
)

print(p)
