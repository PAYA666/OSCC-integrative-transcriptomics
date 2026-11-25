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

ppi_file <- file.path(
  proj_dir,
  "string_interactions (5).tsv"
)

stopifnot(
  file.exists(ppi_file)
)

pkgs <- c(
  "igraph",
  "ggplot2",
  "ggrepel"
)

for (p in pkgs) {

  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }

}

library(igraph)
library(ggplot2)
library(ggrepel)

cand_file <- file.path(
  fig_dir,
  "Figure4B_SourceData.tsv"
)

if (file.exists(cand_file)) {

  cand <- read.delim(
    cand_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

} else {

  cand <- data.frame(

    gene_name = c(
      "SFRP4",
      "PDGFRA",
      "SFRP2",
      "TK1",
      "OGN",
      "MYBL2",
      "FNDC1"
    ),

    Module = c(
      "turquoise",
      "turquoise",
      "turquoise",
      "black",
      "turquoise",
      "black",
      "turquoise"
    ),

    PPI_Degree = c(
      14,
      9,
      10,
      6,
      9,
      5,
      6
    ),

    stringsAsFactors = FALSE
  )
}

seed_genes <- cand$gene_name

ppi_raw <- read.delim(
  ppi_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "STRING table:",
  nrow(ppi_raw),
  "edges x",
  ncol(ppi_raw),
  "columns\n"
)

cat("\nSTRING columns:\n")
print(
  colnames(ppi_raw)
)

node1_candidates <- c(
  "node1_preferred_name",
  "preferredName_A",
  "preferredName.A",
  "node1",
  "#node1"
)

node2_candidates <- c(
  "node2_preferred_name",
  "preferredName_B",
  "preferredName.B",
  "node2"
)

node1_col <- node1_candidates[
  node1_candidates %in%
    colnames(ppi_raw)
][1]

node2_col <- node2_candidates[
  node2_candidates %in%
    colnames(ppi_raw)
][1]

if (
  is.na(node1_col) ||
  is.na(node2_col)
) {

  stop(
    paste0(
      "Could not automatically identify STRING gene-name columns.\n",
      "Please inspect the printed column names."
    )
  )
}

cat(
  "\nUsing node columns:\n",
  node1_col,
  "\n",
  node2_col,
  "\n"
)

score_candidates <- c(
  "score",
  "combined_score",
  "combined.score",
  "combinedScore"
)

score_col <- score_candidates[
  score_candidates %in%
    colnames(ppi_raw)
][1]

edges <- data.frame(

  from = toupper(
    trimws(
      as.character(
        ppi_raw[[node1_col]]
      )
    )
  ),

  to = toupper(
    trimws(
      as.character(
        ppi_raw[[node2_col]]
      )
    )
  ),

  stringsAsFactors = FALSE
)

if (!is.na(score_col)) {

  edges$score <- suppressWarnings(
    as.numeric(
      ppi_raw[[score_col]]
    )
  )

  if (
    max(
      edges$score,
      na.rm = TRUE
    ) > 1.5
  ) {

    edges$score <-
      edges$score / 1000
  }

} else {

  edges$score <- 1
}

edges <- edges[
  !is.na(edges$from) &
  !is.na(edges$to) &
  edges$from != "" &
  edges$to != "" &
  edges$from != edges$to,
]

edges$key <- apply(
  cbind(
    edges$from,
    edges$to
  ),
  1,
  function(x) {
    paste(
      sort(x),
      collapse = "|||"
    )
  }
)

edges <- edges[
  order(
    -edges$score
  ),
]

edges <- edges[
  !duplicated(edges$key),
]

edges$key <- NULL

g_full <- graph_from_data_frame(
  edges[
    ,
    c(
      "from",
      "to",
      "score"
    )
  ],
  directed = FALSE
)

cat(
  "\nFull PPI graph:",
  vcount(g_full),
  "nodes /",
  ecount(g_full),
  "edges\n"
)

seed_present <- intersect(
  toupper(seed_genes),
  V(g_full)$name
)

seed_missing <- setdiff(
  toupper(seed_genes),
  V(g_full)$name
)

cat("\nSeed genes present:\n")
print(seed_present)

if (length(seed_missing) > 0) {

  cat("\nSeed genes missing from STRING network:\n")
  print(seed_missing)

}

neighbor_list <- lapply(
  seed_present,
  function(g) {

    names(
      neighbors(
        g_full,
        g,
        mode = "all"
      )
    )

  }
)

names(neighbor_list) <-
  seed_present

all_neighbors <- unique(
  unlist(
    neighbor_list
  )
)

nonseed_neighbors <- setdiff(
  all_neighbors,
  seed_present
)

neighbor_stats <- lapply(
  nonseed_neighbors,
  function(gene) {

    seed_links <- sum(
      vapply(
        seed_present,
        function(seed) {

          are_adjacent(
            g_full,
            gene,
            seed
          )

        },
        logical(1)
      )
    )

    full_degree <- degree(
      g_full,
      v = gene
    )

    incident_edges <- incident(
      g_full,
      gene
    )

    mean_score <- mean(
      E(g_full)[incident_edges]$score,
      na.rm = TRUE
    )

    data.frame(

      gene_name = gene,

      CandidateLinks =
        seed_links,

      FullDegree =
        as.numeric(
          full_degree
        ),

      MeanScore =
        mean_score,

      stringsAsFactors = FALSE
    )

  }
)

neighbor_stats <- do.call(
  rbind,
  neighbor_stats
)

neighbor_stats <- neighbor_stats[
  order(
    -neighbor_stats$CandidateLinks,
    -neighbor_stats$FullDegree,
    -neighbor_stats$MeanScore
  ),
]

n_neighbors_keep <- 12

top_neighbors <- head(
  neighbor_stats$gene_name,
  n_neighbors_keep
)

local_nodes <- unique(
  c(
    seed_present,
    top_neighbors
  )
)

g_local <- induced_subgraph(
  g_full,
  vids = local_nodes
)

isolated <- names(
  which(
    degree(g_local) == 0
  )
)

if (length(isolated) > 0) {

  g_local <- delete_vertices(
    g_local,
    isolated
  )

}

cat(
  "\nLocal PPI graph:",
  vcount(g_local),
  "nodes /",
  ecount(g_local),
  "edges\n"
)

set.seed(2026)

lay <- layout_with_fr(
  g_local,
  weights = E(g_local)$score,
  niter = 2000
)

node_df <- data.frame(

  gene_name =
    V(g_local)$name,

  x =
    lay[, 1],

  y =
    lay[, 2],

  stringsAsFactors = FALSE
)

node_df$IsCandidate <-
  node_df$gene_name %in%
  seed_present

cand_map <- cand

cand_map$gene_name <-
  toupper(
    cand_map$gene_name
  )

node_df$Module <- cand_map$Module[
  match(
    node_df$gene_name,
    cand_map$gene_name
  )
]

node_df$NodeGroup <- "First-order interactor"

node_df$NodeGroup[
  node_df$IsCandidate &
  node_df$Module == "turquoise"
] <- "Turquoise candidate"

node_df$NodeGroup[
  node_df$IsCandidate &
  node_df$Module == "black"
] <- "Black candidate"

node_df$LocalDegree <- degree(
  g_local,
  v = node_df$gene_name
)

node_df$FullDegree <- degree(
  g_full,
  v = node_df$gene_name
)

bridge_genes <- character(0)

if (nrow(neighbor_stats) > 0) {

  bridge_genes <- neighbor_stats$gene_name[
    neighbor_stats$CandidateLinks >= 2
  ]

}

node_df$Label <- ifelse(

  node_df$IsCandidate |
  node_df$gene_name %in% bridge_genes,

  node_df$gene_name,

  NA_character_

)

edge_df <- as_data_frame(
  g_local,
  what = "edges"
)

edge_df$x <-
  node_df$x[
    match(
      edge_df$from,
      node_df$gene_name
    )
  ]

edge_df$y <-
  node_df$y[
    match(
      edge_df$from,
      node_df$gene_name
    )
  ]

edge_df$xend <-
  node_df$x[
    match(
      edge_df$to,
      node_df$gene_name
    )
  ]

edge_df$yend <-
  node_df$y[
    match(
      edge_df$to,
      node_df$gene_name
    )
  ]

if (
  diff(
    range(
      edge_df$score,
      na.rm = TRUE
    )
  ) == 0
) {

  edge_df$edge_width <- 0.7
  edge_df$edge_alpha <- 0.45

} else {

  edge_df$edge_width <- 0.35 +
    1.15 *
    (
      edge_df$score -
      min(
        edge_df$score,
        na.rm = TRUE
      )
    ) /
    diff(
      range(
        edge_df$score,
        na.rm = TRUE
      )
    )

  edge_df$edge_alpha <- 0.25 +
    0.50 *
    (
      edge_df$score -
      min(
        edge_df$score,
        na.rm = TRUE
      )
    ) /
    diff(
      range(
        edge_df$score,
        na.rm = TRUE
      )
    )

}

node_colors <- c(

  "Turquoise candidate" =
    "#4C8D8A",

  "Black candidate" =
    "#444444",

  "First-order interactor" =
    "#D3D3D3"
)

p4c <- ggplot() +

  geom_segment(

    data = edge_df,

    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),

    linewidth =
      edge_df$edge_width,

    alpha =
      edge_df$edge_alpha,

    color = "#A7A7A7",

    lineend = "round"

  ) +

  geom_point(

    data = node_df[
      !node_df$IsCandidate,
    ],

    aes(
      x = x,
      y = y
    ),

    shape = 21,

    size = 4.0,

    fill = "#D3D3D3",

    color = "#777777",

    stroke = 0.45

  ) +

  geom_point(

    data = node_df[
      node_df$IsCandidate,
    ],

    aes(
      x = x,
      y = y,
      fill = NodeGroup,
      size = FullDegree
    ),

    shape = 21,

    color = "black",

    stroke = 0.85

  ) +

  geom_text_repel(

    data = node_df[
      !is.na(node_df$Label),
    ],

    aes(
      x = x,
      y = y,
      label = Label,
      fontface = ifelse(
        IsCandidate,
        "bold.italic",
        "italic"
      )
    ),

    size = 3.3,

    color = "black",

    box.padding = 0.5,

    point.padding = 0.35,

    min.segment.length = 0,

    segment.color = "#777777",

    segment.size = 0.3,

    max.overlaps = Inf,

    seed = 2026

  ) +

  scale_fill_manual(

    values = node_colors[
      c(
        "Turquoise candidate",
        "Black candidate"
      )
    ],

    name = NULL

  ) +

  scale_size_continuous(

    range = c(
      6.0,
      10.0
    ),

    name = "PPI degree"

  ) +

  coord_equal(
    clip = "off"
  ) +

  theme_void(
    base_size = 11.5,
    base_family = "sans"
  ) +

  theme(

    legend.position = "right",

    legend.title = element_text(
      size = 9.5
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      25,
      35,
      25,
      35
    )

  )

print(
  p4c
)

ggsave(

  filename = file.path(
    fig_dir,
    "Figure4C_Local_STRING_PPI_PUBLICATION.pdf"
  ),

  plot = p4c,

  width = 8.0,

  height = 7.0,

  units = "in",

  device = cairo_pdf

)

ggsave(

  filename = file.path(
    fig_dir,
    "Figure4C_Local_STRING_PPI_PUBLICATION.tiff"
  ),

  plot = p4c,

  width = 8.0,

  height = 7.0,

  units = "in",

  dpi = 600,

  compression = "lzw"

)

ggsave(

  filename = file.path(
    fig_dir,
    "Figure4C_Local_STRING_PPI_PUBLICATION.png"
  ),

  plot = p4c,

  width = 8.0,

  height = 7.0,

  units = "in",

  dpi = 200

)

write.table(

  node_df,

  file.path(
    fig_dir,
    "Figure4C_SourceData_Nodes.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

write.table(

  edge_df,

  file.path(
    fig_dir,
    "Figure4C_SourceData_Edges.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE

)

cat("\n========================================\n")
cat("Figure 4C local STRING-PPI finished\n")
cat("========================================\n")

cat(
  "Full network:",
  vcount(g_full),
  "nodes /",
  ecount(g_full),
  "edges\n"
)

cat(
  "Local network:",
  vcount(g_local),
  "nodes /",
  ecount(g_local),
  "edges\n"
)

cat("\nCandidate genes in local graph:\n")

print(
  node_df[
    node_df$IsCandidate,
    c(
      "gene_name",
      "Module",
      "FullDegree",
      "LocalDegree"
    )
  ]
)

cat("\nSelected first-order interactors:\n")

print(
  neighbor_stats[
    neighbor_stats$gene_name %in%
      top_neighbors,
    ]
)
