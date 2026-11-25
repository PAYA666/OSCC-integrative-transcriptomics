rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages(
    "ggplot2",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages(
    "ggrepel",
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

part3_dir <- file.path(project_root, "Part3_Cisplatin")

step1_file <- file.path(
  part3_dir,
  "Step1_FINAL_V4_OFFLINE",
  "Part3_Step1_FINAL_V4_AllResults.rds"
)

out_dir <- file.path(
  part3_dir,
  "Figures_FINAL",
  "Figure1"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(step1_file)) {
  stop(
    paste0(
      "Cannot find Step1 RDS:\n",
      step1_file
    )
  )
}

sample_ids <- c(
  "GSM5920580",
  "GSM5920581",
  "GSM5920582",
  "GSM5920583",
  "GSM5920584",
  "GSM5920585"
)

group <- c(
  "Resistant",
  "Resistant",
  "Resistant",
  "Parental",
  "Parental",
  "Parental"
)

group <- factor(
  group,
  levels = c(
    "Parental",
    "Resistant"
  )
)

step1 <- readRDS(
  step1_file
)

resolve_sample_columns <- function(column_names, sample_ids) {

  if (is.null(column_names)) {
    return(NULL)
  }

  if (all(sample_ids %in% column_names)) {
    return(sample_ids)
  }

  resolved <- character(
    length(sample_ids)
  )

  for (i in seq_along(sample_ids)) {

    hit <- grep(
      sample_ids[i],
      column_names,
      fixed = TRUE
    )

    if (length(hit) != 1) {
      return(NULL)
    }

    resolved[i] <- column_names[hit]
  }

  return(resolved)
}

find_expression_matrices <- function(x, sample_ids) {

  results <- list()

  walk_object <- function(obj, object_path) {

    if (is.matrix(obj) || is.data.frame(obj)) {

      if (nrow(obj) > 1000 && ncol(obj) >= 6) {

        sample_columns <- resolve_sample_columns(
          colnames(obj),
          sample_ids
        )

        if (!is.null(sample_columns)) {

          results[[length(results) + 1L]] <<- list(
            path = object_path,
            nrow = nrow(obj),
            ncol = ncol(obj),
            sample_columns = sample_columns,
            object = obj
          )
        }
      }
    }

    if (is.list(obj) && !is.data.frame(obj)) {

      object_names <- names(obj)

      if (is.null(object_names)) {
        object_names <- rep(
          "",
          length(obj)
        )
      }

      for (i in seq_along(obj)) {

        current_name <- object_names[i]

        if (
          is.na(current_name) ||
          current_name == ""
        ) {

          child_path <- paste0(
            object_path,
            "[[",
            i,
            "]]"
          )

        } else {

          child_path <- paste0(
            object_path,
            "$",
            current_name
          )
        }

        walk_object(
          obj[[i]],
          child_path
        )
      }
    }
  }

  walk_object(
    step1,
    "step1"
  )

  return(results)
}

expression_candidates <- find_expression_matrices(
  step1,
  sample_ids
)

if (length(expression_candidates) == 0) {

  stop(
    paste0(
      "No expression matrix containing all six GSE197561 samples ",
      "was found in the Step1 RDS."
    )
  )
}

candidate_summary <- data.frame(
  Index = seq_along(expression_candidates),
  Path = vapply(
    expression_candidates,
    function(x) x$path,
    character(1)
  ),
  Rows = vapply(
    expression_candidates,
    function(x) x$nrow,
    numeric(1)
  ),
  Columns = vapply(
    expression_candidates,
    function(x) x$ncol,
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

candidate_summary$Distance_from_26972 <- abs(
  candidate_summary$Rows - 26972
)

best_index <- candidate_summary$Index[
  which.min(
    candidate_summary$Distance_from_26972
  )
]

selected <- expression_candidates[[best_index]]

expr <- selected$object[
  ,
  selected$sample_columns,
  drop = FALSE
]

expr <- as.matrix(expr)

suppressWarnings(
  storage.mode(expr) <- "numeric"
)

colnames(expr) <- sample_ids

finite_rows <- apply(
  expr,
  1,
  function(x) {
    all(is.finite(x))
  }
)

expr <- expr[
  finite_rows,
  ,
  drop = FALSE
]

row_sd <- apply(
  expr,
  1,
  sd
)

expr <- expr[
  is.finite(row_sd) &
    row_sd > 0,
  ,
  drop = FALSE
]

if (nrow(expr) < 1000) {
  stop("Too few valid genes remained for PCA.")
}

pca <- prcomp(
  t(expr),
  center = TRUE,
  scale. = FALSE
)

variance_percent <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100

pca_df <- data.frame(
  Sample = sample_ids,
  Group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_df$Label <- c(
  "R1",
  "R2",
  "R3",
  "P1",
  "P2",
  "P3"
)

group_colors <- c(
  "Parental" = "#3C6E9F",
  "Resistant" = "#B64D43"
)

p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = Group
  )
) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "#E5E5E5"
  ) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "#E5E5E5"
  ) +

  geom_point(
    size = 3.8,
    alpha = 0.95
  ) +

  ggrepel::geom_text_repel(
    aes(
      label = Label
    ),
    size = 3.2,
    colour = "#222222",
    box.padding = 0.35,
    point.padding = 0.30,
    segment.size = 0.30,
    segment.colour = "#777777",
    min.segment.length = 0,
    seed = 20260815,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = group_colors
  ) +

  labs(
    title = "GSE197561",
    subtitle = "Acquired cisplatin resistance",
    x = sprintf(
      "PC1 (%.1f%%)",
      variance_percent[1]
    ),
    y = sprintf(
      "PC2 (%.1f%%)",
      variance_percent[2]
    ),
    colour = NULL
  ) +

  theme_classic(
    base_size = 11,
    base_family = "sans"
  ) +

  theme(
    axis.title = element_text(
      size = 10,
      colour = "black"
    ),

    axis.text = element_text(
      size = 9,
      colour = "black"
    ),

    axis.line = element_line(
      linewidth = 0.55,
      colour = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.45,
      colour = "black"
    ),

    plot.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),

    plot.subtitle = element_text(
      size = 9.5,
      colour = "#555555"
    ),

    legend.position = "top",

    legend.justification = "left",

    legend.text = element_text(
      size = 9
    ),

    legend.margin = margin(
      b = -3
    ),

    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

pdf_file <- file.path(
  out_dir,
  "Figure1_Plot1_GSE197561_PCA.pdf"
)

png_file <- file.path(
  out_dir,
  "Figure1_Plot1_GSE197561_PCA.png"
)

tiff_file <- file.path(
  out_dir,
  "Figure1_Plot1_GSE197561_PCA_600dpi.tiff"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 88,
  height = 78,
  units = "mm",
  device = "pdf",
  bg = "white",
  useDingbats = FALSE
)

ggsave(
  filename = png_file,
  plot = p,
  width = 88,
  height = 78,
  units = "mm",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = p,
  width = 88,
  height = 78,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write.table(
  pca_df,
  file = file.path(
    out_dir,
    "Figure1_Plot1_GSE197561_PCA_coordinates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pca_metrics <- data.frame(
  Dataset = "GSE197561",
  Expression_Rows_Used = nrow(expr),
  Samples = ncol(expr),
  PC1_Percent = variance_percent[1],
  PC2_Percent = variance_percent[2],
  Selected_Object_Path = selected$path,
  stringsAsFactors = FALSE
)

write.table(
  pca_metrics,
  file = file.path(
    out_dir,
    "Figure1_Plot1_GSE197561_PCA_metrics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  list(
    Plot = p,
    PCA_Data = pca_df,
    Metrics = pca_metrics,
    Candidate_Matrices = candidate_summary
  ),
  file = file.path(
    out_dir,
    "Figure1_Plot1_GSE197561_PCA_object.rds"
  )
)

print(p)

cat("\n")
cat("============================================\n")
cat("GSE197561 PCA COMPLETE\n")
cat("============================================\n")

print(pca_metrics)

cat(
  "\nPNG:\n",
  png_file,
  "\n"
)
