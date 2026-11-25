project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root
fig_dir  <- file.path(
  proj_dir,
  "Figures_Final"
)

p6_final <- p6 +

  theme(

    plot.background = element_rect(
      fill = "white",
      color = NA
    ),

    panel.background = element_rect(
      fill = "white",
      color = NA
    )

  )

print(
  p6_final
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure6_WorkingModel_PUBLICATION_FINAL.pdf"
  ),
  plot = p6_final,
  width = 11.2,
  height = 7.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure6_WorkingModel_PUBLICATION_FINAL.tiff"
  ),
  plot = p6_final,
  width = 11.2,
  height = 7.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure6_WorkingModel_PUBLICATION_FINAL.png"
  ),
  plot = p6_final,
  width = 11.2,
  height = 7.6,
  units = "in",
  dpi = 250,
  bg = "white"
)

cat(
  "Figure6 white-background version saved.\n"
)
