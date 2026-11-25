project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure3")

score_df <- step4$scores

cor_test <- cor.test(
  score_df$Priority52_Score,
  score_df$Ultra7_Score,
  method = "spearman"
)

rho_text <- paste0(
  "Spearman rho = ",
  round(
    cor_test$estimate,
    2
  ),
  "\nP = ",
  format.pval(
    cor_test$p.value,
    digits = 2
  )
)

p <- ggplot(
  score_df,
  aes(
    x = Priority52_Score,
    y = Ultra7_Score
  )
)+

geom_point(
  size = 2,
  alpha = 0.6
)+

geom_smooth(
  method = "lm",
  se = FALSE
)+

annotate(
  "text",
  x = min(score_df$Priority52_Score),
  y = max(score_df$Ultra7_Score),
  label = rho_text,
  hjust = 0,
  vjust = 1,
  size = 5
)+

theme_classic(
  base_size = 16
)+

labs(
  title =
    "Ultra7 score coherence with the broader resistance program",

  subtitle =
    "TCGA HNSCC cohort",

  x =
    "Priority52 resistance program score",

  y =
    "Ultra7 score"
)

ggsave(
  file.path(
    out_dir,
    "Figure3D_Ultra7_TCGA_score_coherence.png"
  ),

  p,

  width = 85,
  height = 75,
  units = "mm",
  dpi = 600,
  bg = "white"
)

print(p)
