project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

library(igraph)

string_file <- file.path(
  proj_dir,
  "string_interactions (5).tsv"
)

stopifnot(file.exists(string_file))

string_raw <- read.delim(
  string_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("STRING dimensions:\n")
print(dim(string_raw))

cat("\nSTRING columns:\n")
print(colnames(string_raw))

possible_A <- c(
  "#node1",
  "node1",
  "preferredName_A",
  "protein1"
)

possible_B <- c(
  "node2",
  "preferredName_B",
  "protein2"
)

col_A <- intersect(
  possible_A,
  colnames(string_raw)
)

col_B <- intersect(
  possible_B,
  colnames(string_raw)
)

if (length(col_A) == 0 | length(col_B) == 0) {
  stop(
    "Cannot automatically identify STRING node columns. Check colnames(string_raw)."
  )
}

col_A <- col_A[1]
col_B <- col_B[1]

cat("\nUsing node columns:\n")
cat("A =", col_A, "\n")
cat("B =", col_B, "\n")

edges <- data.frame(
  from = string_raw[[col_A]],
  to   = string_raw[[col_B]],
  stringsAsFactors = FALSE
)

edges <- edges[
  !is.na(edges$from) &
  !is.na(edges$to) &
  edges$from != "" &
  edges$to != "" &
  edges$from != edges$to,
]

edge_key <- apply(
  edges,
  1,
  function(x) {
    paste(
      sort(x),
      collapse = "__"
    )
  }
)

edges <- edges[
  !duplicated(edge_key),
]

cat(
  "\nUnique undirected PPI edges:",
  nrow(edges),
  "\n"
)

g <- graph_from_data_frame(
  edges,
  directed = FALSE
)

g <- simplify(
  g,
  remove.multiple = TRUE,
  remove.loops = TRUE
)

cat(
  "PPI nodes:",
  vcount(g),
  "\n"
)

cat(
  "PPI edges:",
  ecount(g),
  "\n"
)

comp <- components(g)

component_table <- sort(
  table(comp$membership),
  decreasing = TRUE
)

cat("\nConnected component sizes:\n")
print(component_table)

largest_component_id <- as.numeric(
  names(component_table)[1]
)

main_nodes <- names(
  comp$membership[
    comp$membership == largest_component_id
  ]
)

g_main <- induced_subgraph(
  g,
  vids = main_nodes
)

cat(
  "\nMain component nodes:",
  vcount(g_main),
  "\n"
)

cat(
  "Main component edges:",
  ecount(g_main),
  "\n"
)

Degree <- degree(
  g,
  mode = "all"
)

Betweenness <- betweenness(
  g,
  directed = FALSE,
  normalized = TRUE
)

Closeness <- closeness(
  g,
  normalized = TRUE
)

ppi_topology <- data.frame(
  gene_name = names(Degree),

  Degree = as.numeric(Degree),

  Betweenness = as.numeric(
    Betweenness[names(Degree)]
  ),

  Closeness = as.numeric(
    Closeness[names(Degree)]
  ),

  stringsAsFactors = FALSE
)

ppi_topology$Component <- comp$membership[
  ppi_topology$gene_name
]

ppi_topology$MainComponent <-
  ppi_topology$gene_name %in%
  main_nodes

ppi_topology$DegreeRank <- rank(
  -ppi_topology$Degree,
  ties.method = "min"
)

ppi_topology$BetweennessRank <- rank(
  -ppi_topology$Betweenness,
  ties.method = "min"
)

ppi_topology$ClosenessRank <- rank(
  -ppi_topology$Closeness,
  ties.method = "min"
)

ppi_topology$MeanTopologyRank <- rowMeans(
  ppi_topology[
    ,
    c(
      "DegreeRank",
      "BetweennessRank",
      "ClosenessRank"
    )
  ],
  na.rm = TRUE
)

ppi_topology <- ppi_topology[
  order(
    ppi_topology$MeanTopologyRank,
    ppi_topology$DegreeRank
  ),
]

wgcna_file <- file.path(
  proj_dir,
  "Part2_Smoking_WGCNA_GeneLevel_MM_GS.rds"
)

stopifnot(file.exists(wgcna_file))

wgcna_gene <- readRDS(
  wgcna_file
)

integrated <- merge(
  ppi_topology,
  wgcna_gene,
  by = "gene_name",
  all.x = TRUE
)

integrated$WGCNA_Core <- (
  !is.na(integrated$MM) &
  abs(integrated$MM) >= 0.70 &
  integrated$Smoking_DEG &
  integrated$module_color %in%
    c(
      "turquoise",
      "black",
      "red"
    )
)

integrated$PPI_Top20 <- (
  integrated$MeanTopologyRank <= 20
)

integrated$MultiEvidence <- (
  integrated$WGCNA_Core &
  integrated$PPI_Top20
)

integrated$EvidenceCount <- 0L

integrated$EvidenceCount <-
  integrated$EvidenceCount +
  ifelse(
    !is.na(integrated$Smoking_DEG) &
    integrated$Smoking_DEG,
    1L,
    0L
  )

integrated$EvidenceCount <-
  integrated$EvidenceCount +
  ifelse(
    !is.na(integrated$module_color) &
    integrated$module_color %in%
      c(
        "turquoise",
        "black",
        "red"
      ),
    1L,
    0L
  )

integrated$EvidenceCount <-
  integrated$EvidenceCount +
  ifelse(
    !is.na(integrated$MM) &
    abs(integrated$MM) >= 0.70,
    1L,
    0L
  )

integrated$EvidenceCount <-
  integrated$EvidenceCount +
  ifelse(
    integrated$PPI_Top20,
    1L,
    0L
  )

integrated <- integrated[
  order(
    -integrated$EvidenceCount,
    integrated$MeanTopologyRank,
    -abs(integrated$MM)
  ),
]

ppi_top20 <- integrated[
  integrated$PPI_Top20,
]

ppi_top20 <- ppi_top20[
  order(
    ppi_top20$MeanTopologyRank
  ),
]

multi_candidates <- integrated[
  integrated$MultiEvidence,
]

multi_candidates <- multi_candidates[
  order(
    multi_candidates$MeanTopologyRank
  ),
]

write.table(
  ppi_topology,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_Topology.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  integrated,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_WGCNA_Integrated.rds"
  )
)

write.table(
  integrated,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_WGCNA_Integrated.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ppi_top20,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_Top20.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  multi_candidates,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_WGCNA_MultiEvidence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 10B finished\n")
cat("========================================\n")

cat(
  "PPI nodes:",
  nrow(ppi_topology),
  "\n"
)

cat(
  "Main component:",
  sum(ppi_topology$MainComponent),
  "\n"
)

cat(
  "PPI Top20:",
  nrow(ppi_top20),
  "\n"
)

cat(
  "WGCNA-core + PPI-Top20:",
  nrow(multi_candidates),
  "\n"
)
