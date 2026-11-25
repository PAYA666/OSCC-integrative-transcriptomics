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

deg_key <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_DESeq2_FINAL.rds"
  )
)

status_tab <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_StatusPattern_FINAL.rds"
  )
)

fib_tab <- read.delim(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_FibroblastCorrelation.tsv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

network_candidates <- c(
  "SFRP4",
  "PDGFRA",
  "SFRP2",
  "TK1",
  "OGN",
  "MYBL2",
  "FNDC1"
)

gene_order <- c(

  "CYP1A1",
  "CES1",

  "SFRP2",
  "PDGFRA",
  "SFRP4",
  "FNDC1",
  "OGN",

  "TK1",
  "MYBL2"
)

gene_info <- data.frame(

  gene_name = gene_order,

  Program = c(
    "Exposure-response",
    "Exposure-response",

    "Stromal-associated",
    "Stromal-associated",
    "Stromal-associated",
    "Stromal-associated",
    "Stromal-associated",

    "Proliferative",
    "Proliferative"
  ),

  stringsAsFactors = FALSE
)

integrated <- gene_info

integrated$Smoking_DEG <- deg_key$padj[
  match(
    integrated$gene_name,
    deg_key$gene_name
  )
] < 0.05

integrated$DEG_log2FC <- deg_key$log2FoldChange[
  match(
    integrated$gene_name,
    deg_key$gene_name
  )
]

integrated$DEG_FDR <- deg_key$padj[
  match(
    integrated$gene_name,
    deg_key$gene_name
  )
]

integrated$Current_vs_Never <-
  status_tab$Current_vs_Never_FDR[
    match(
      integrated$gene_name,
      status_tab$gene_name
    )
  ] < 0.05

integrated$Former_vs_Never <-
  status_tab$Former_vs_Never_FDR[
    match(
      integrated$gene_name,
      status_tab$gene_name
    )
  ] < 0.05

integrated$PackYears_Rho <-
  status_tab$PackYears_Rho[
    match(
      integrated$gene_name,
      status_tab$gene_name
    )
  ]

integrated$PackYears_FDR <-
  status_tab$PackYears_FDR[
    match(
      integrated$gene_name,
      status_tab$gene_name
    )
  ]

integrated$Fibroblast_Rho <-
  fib_tab$Spearman_Rho[
    match(
      integrated$gene_name,
      fib_tab$gene_name
    )
  ]

integrated$Fibroblast_FDR <-
  fib_tab$FDR[
    match(
      integrated$gene_name,
      fib_tab$gene_name
    )
  ]

integrated$Fibroblast_Associated <-
  !is.na(integrated$Fibroblast_Rho) &
  integrated$Fibroblast_Rho > 0.50 &
  integrated$Fibroblast_FDR < 0.05

integrated$WGCNA_Core <-
  integrated$gene_name %in%
  network_candidates

integrated$PPI_Top20 <-
  integrated$gene_name %in%
  network_candidates

integrated$Multi_Evidence <-
  integrated$WGCNA_Core &
  integrated$PPI_Top20

print(
  integrated
)

evidence_columns <- c(
  "Smoking_DEG",
  "Current_vs_Never",
  "Former_vs_Never",
  "Fibroblast_Associated",
  "WGCNA_Core",
  "PPI_Top20",
  "Multi_Evidence"
)

evidence_labels <- c(

  "Smoking DEG",

  "Current vs Never",

  "Former vs Never",

  "Fibroblast\nassociated",

  "WGCNA core",

  "PPI Top20",

  "Multi-evidence"

)

matrix_list <- lapply(
  seq_len(nrow(integrated)),
  function(i) {

    data.frame(

      gene_name =
        integrated$gene_name[i],

      Program =
        integrated$Program[i],

      Evidence =
        evidence_columns,

      Present =
        as.logical(
          integrated[
            i,
            evidence_columns
          ]
        ),

      stringsAsFactors = FALSE
    )
  }
)

matrix_dat <- do.call(
  rbind,
  matrix_list
)

matrix_dat$EvidenceType <- "Network"

matrix_dat$EvidenceType[
  matrix_dat$Evidence %in%
    c(
      "Smoking_DEG",
      "Current_vs_Never",
      "Former_vs_Never"
    )
] <- "Smoking association"

matrix_dat$EvidenceType[
  matrix_dat$Evidence ==
    "Fibroblast_Associated"
] <- "Stromal context"

matrix_dat$EvidenceType[
  matrix_dat$Evidence ==
    "Multi_Evidence"
] <- "Integrated evidence"

matrix_dat$gene_name <- factor(
  matrix_dat$gene_name,
  levels = rev(
    gene_order
  )
)

matrix_dat$Evidence <- factor(
  matrix_dat$Evidence,
  levels = evidence_columns,
  labels = evidence_labels
)

matrix_dat$EvidenceType <- factor(
  matrix_dat$EvidenceType,
  levels = c(
    "Smoking association",
    "Stromal context",
    "Network",
    "Integrated evidence"
  )
)

evidence_colors <- c(

  "Smoking association" =
    "#B75B52",

  "Stromal context" =
    "#4E8077",

  "Network" =
    "#596B8C",

  "Integrated evidence" =
    "#7A557D"
)

p5a <- ggplot(
  matrix_dat,
  aes(
    x = Evidence,
    y = gene_name
  )
) +

  geom_point(
    shape = 21,
    size = 5.6,
    fill = "white",
    color = "#D2D2D2",
    stroke = 0.7
  ) +

  geom_point(
    data = matrix_dat[
      matrix_dat$Present,
    ],
    aes(
      fill = EvidenceType
    ),
    shape = 21,
    size = 5.6,
    color = "black",
    stroke = 0.55
  ) +

  geom_vline(
    xintercept = c(
      3.5,
      4.5,
      6.5
    ),
    color = "#D8D8D8",
    linewidth = 0.45
  ) +

  scale_fill_manual(
    values = evidence_colors,
    name = NULL
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(
    base_size = 11.5,
    base_family = "sans"
  ) +

  theme(

    axis.line =
      element_blank(),

    axis.ticks =
      element_blank(),

    axis.text.y =
      element_text(
        size = 10.5,
        face = "italic",
        color = "black"
      ),

    axis.text.x =
      element_text(
        size = 9.6,
        color = "black",
        angle = 35,
        hjust = 1,
        vjust = 1
      ),

    legend.position =
      "top",

    legend.direction =
      "horizontal",

    legend.text =
      element_text(
        size = 9.2
      ),

    panel.grid =
      element_blank(),

    plot.margin =
      margin(
        15,
        15,
        15,
        12
      )
  )

program_annotation <- data.frame(

  gene_name = factor(
    c(
      "CES1",
      "OGN",
      "MYBL2"
    ),
    levels = rev(gene_order)
  ),

  ProgramLabel = c(
    "Exposure-response",
    "Stromal-associated",
    "Proliferative"
  ),

  stringsAsFactors = FALSE
)

p5a <- p5a +

  geom_text(
    data = program_annotation,
    aes(
      x = 7.75,
      y = gene_name,
      label = ProgramLabel
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.3,
    fontface = "bold",
    color = "#555555"
  ) +

  coord_cartesian(
    xlim = c(
      0.6,
      9.1
    ),
    clip = "off"
  )

print(
  p5a
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure5A_EvidenceConvergenceMatrix_PUBLICATION.pdf"
  ),
  plot = p5a,
  width = 9.2,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure5A_EvidenceConvergenceMatrix_PUBLICATION.tiff"
  ),
  plot = p5a,
  width = 9.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    fig_dir,
    "Figure5A_EvidenceConvergenceMatrix_PUBLICATION.png"
  ),
  plot = p5a,
  width = 9.2,
  height = 5.8,
  units = "in",
  dpi = 200
)

write.table(
  integrated,
  file.path(
    fig_dir,
    "Figure5A_SourceData_IntegratedGenes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  matrix_dat,
  file.path(
    fig_dir,
    "Figure5A_SourceData_EvidenceMatrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Figure 5A finished\n")
cat("========================================\n")

print(
  integrated[
    ,
    c(
      "gene_name",
      "Program",
      "Smoking_DEG",
      "Current_vs_Never",
      "Former_vs_Never",
      "Fibroblast_Associated",
      "WGCNA_Core",
      "PPI_Top20",
      "Multi_Evidence"
    )
  ]
)
