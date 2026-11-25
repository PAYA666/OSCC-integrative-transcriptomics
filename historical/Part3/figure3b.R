project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(reshape2)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure3")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cor_mat <- step4$score_correlations

print(cor_mat)

plot_df <- melt(
  cor_mat
)

colnames(plot_df) <- c(
  "Program1",
  "Program2",
  "Correlation"
)

p <- ggplot(
  plot_df,
  aes(
    x=Program1,
    y=Program2,
    fill=Correlation
  )
)+

geom_tile(
  color="white"
)+

geom_text(
  aes(
    label=
      sprintf(
        "%.2f",
        Correlation
      )
  ),
  size=5
)+

scale_fill_gradient2(
  low="#4575B4",
  mid="white",
  high="#B84D43",
  midpoint=0,
  limits=c(-1,1)
)+

theme_classic(
  base_size=13
)+

theme(
  axis.text.x=
    element_text(
      angle=45,
      hjust=1
    )
)+

labs(
  title=
    "Coherence among frozen cisplatin-resistance programs",

  x=NULL,
  y=NULL,
  fill="Spearman\ncorrelation"
)

ggsave(
  file.path(
    out_dir,
    "Figure3B_Resistance_program_correlation.png"
  ),

  p,

  width=85,
  height=75,
  units="mm",
  dpi=600,
  bg="white"
)

print(p)
