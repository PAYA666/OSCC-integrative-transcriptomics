project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

if (!requireNamespace("STRINGdb", quietly = TRUE)) {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  BiocManager::install(
    "STRINGdb",
    ask = FALSE,
    update = FALSE
  )
}

library(STRINGdb)

target_file <- file.path(
  proj_dir,
  "Part2_Smoking_WGCNA_TargetModules_DEGs.tsv"
)

stopifnot(file.exists(target_file))

target_deg <- read.delim(
  target_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

core_file <- file.path(
  proj_dir,
  "Part2_Smoking_WGCNA_CoreCandidateGenes.rds"
)

stopifnot(file.exists(core_file))

candidate_core <- readRDS(core_file)

count_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_raw_count_matrix.rds"
)

stopifnot(file.exists(count_file))

counts_raw <- readRDS(count_file)

gene_annotation <- unique(
  counts_raw[, c(
    "gene_id",
    "gene_name",
    "gene_type"
  )]
)

idx <- match(
  target_deg$gene_id,
  gene_annotation$gene_id
)

target_deg$gene_type <-
  gene_annotation$gene_type[idx]

ppi_input <- target_deg[
  !is.na(target_deg$gene_name) &
  target_deg$gene_name != "" &
  target_deg$gene_type == "protein_coding",
]

ppi_input <- ppi_input[
  order(
    ppi_input$DEG_padj,
    -abs(ppi_input$log2FoldChange)
  ),
]

ppi_input <- ppi_input[
  !duplicated(ppi_input$gene_name),
]

write.table(
  ppi_input,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_Input_TargetModule_DEGs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ppi_input$gene_name,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_Input_GeneSymbols.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

string_db <- STRINGdb$new(
  species = 9606,
  score_threshold = 400
)

string_input <- data.frame(
  gene = ppi_input$gene_name,
  stringsAsFactors = FALSE
)

string_mapped <- string_db$map(
  string_input,
  "gene",
  removeUnmappedRows = FALSE
)

mapped_ok <- (
  !is.na(string_mapped$STRING_id) &
  string_mapped$STRING_id != ""
)

mapped_genes <- string_mapped[
  mapped_ok,
]

unmapped_genes <- string_mapped[
  !mapped_ok,
]

ppi_edges_raw <- string_db$get_interactions(
  mapped_genes$STRING_id
)

candidate_string_ids <- unique(
  mapped_genes$STRING_id
)

ppi_edges <- ppi_edges_raw[
  ppi_edges_raw$from %in% candidate_string_ids &
  ppi_edges_raw$to %in% candidate_string_ids,
]

id_map <- mapped_genes[
  ,
  c(
    "STRING_id",
    "gene"
  )
]

id_map <- id_map[
  !duplicated(id_map$STRING_id),
]

ppi_edges$from_gene <-
  id_map$gene[
    match(
      ppi_edges$from,
      id_map$STRING_id
    )
  ]

ppi_edges$to_gene <-
  id_map$gene[
    match(
      ppi_edges$to,
      id_map$STRING_id
    )
  ]

ppi_edges <- ppi_edges[
  !is.na(ppi_edges$from_gene) &
  !is.na(ppi_edges$to_gene) &
  ppi_edges$from_gene != ppi_edges$to_gene,
]

edge_key <- apply(
  ppi_edges[, c(
    "from_gene",
    "to_gene"
  )],
  1,
  function(x) {
    paste(
      sort(x),
      collapse = "__"
    )
  }
)

ppi_edges <- ppi_edges[
  !duplicated(edge_key),
]

write.table(
  mapped_genes,
  file.path(
    proj_dir,
    "Part2_Smoking_STRING_MappedGenes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  unmapped_genes,
  file.path(
    proj_dir,
    "Part2_Smoking_STRING_UnmappedGenes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ppi_edges,
  file.path(
    proj_dir,
    "Part2_Smoking_STRING_PPI_Edges.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cytoscape_edges <- ppi_edges[
  ,
  c(
    "from_gene",
    "to_gene",
    "combined_score"
  )
]

colnames(cytoscape_edges) <- c(
  "Source",
  "Target",
  "STRING_score"
)

write.table(
  cytoscape_edges,
  file.path(
    proj_dir,
    "Part2_Smoking_STRING_Cytoscape_Edges.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  ppi_input,
  file.path(
    proj_dir,
    "Part2_Smoking_PPI_Input_TargetModule_DEGs.rds"
  )
)

cat("\n========================================\n")
cat("Step 10A STRING PPI finished\n")
cat("========================================\n")

cat(
  "Target-module Smoking-DEGs:",
  nrow(target_deg),
  "\n"
)

cat(
  "Protein-coding PPI input genes:",
  nrow(ppi_input),
  "\n"
)

cat(
  "Mapped to STRING:",
  nrow(mapped_genes),
  "\n"
)

cat(
  "Unmapped:",
  nrow(unmapped_genes),
  "\n"
)

cat(
  "PPI edges:",
  nrow(ppi_edges),
  "\n"
)
