project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

module_file <- file.path(
  proj_dir,
  "TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv"
)

stopifnot(file.exists(module_file))

gene_module <- read.delim(
  module_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Gene-module table:",
  dim(gene_module),
  "\n"
)

deg_file <- file.path(
  proj_dir,
  "Part2_Smoking_DEG_Annotated.rds"
)

stopifnot(file.exists(deg_file))

deg <- readRDS(deg_file)

required_module_cols <- c(
  "gene_id",
  "module_numeric",
  "module_color",
  "gene_name"
)

stopifnot(
  all(required_module_cols %in% colnames(gene_module))
)

idx <- match(
  gene_module$gene_id,
  deg$gene_id
)

gene_module$log2FoldChange <-
  deg$log2FoldChange[idx]

gene_module$padj <-
  deg$padj[idx]

gene_module$DEG_status <-
  deg$DEG_status[idx]

gene_module$Smoking_DEG <- (
  !is.na(gene_module$padj) &
  gene_module$padj < 0.05
)

gene_module$Smoking_Strong_DEG <- (
  !is.na(gene_module$padj) &
  gene_module$padj < 0.05 &
  abs(gene_module$log2FoldChange) >= 1
)

cat(
  "WGCNA genes matched to DEG results:",
  sum(!is.na(gene_module$padj)),
  "/",
  nrow(gene_module),
  "\n"
)

cat(
  "Smoking-DEGs inside WGCNA network:",
  sum(gene_module$Smoking_DEG),
  "\n"
)

cat(
  "Strong Smoking-DEGs inside WGCNA network:",
  sum(gene_module$Smoking_Strong_DEG),
  "\n"
)

modules <- sort(
  unique(gene_module$module_numeric)
)

enrichment_list <- lapply(
  modules,
  function(mod) {

    in_module <- (
      gene_module$module_numeric == mod
    )

    is_deg <- gene_module$Smoking_DEG

    a <- sum(in_module & is_deg)
    b <- sum(in_module & !is_deg)

    c <- sum(!in_module & is_deg)
    d <- sum(!in_module & !is_deg)

    fisher_mat <- matrix(
      c(a, b, c, d),
      nrow = 2,
      byrow = TRUE
    )

    ft <- fisher.test(
      fisher_mat,
      alternative = "greater"
    )

    module_color_value <- unique(
      gene_module$module_color[in_module]
    )

    module_color_value <- paste(
      module_color_value,
      collapse = ";"
    )

    data.frame(
      module_numeric = mod,
      module_color = module_color_value,
      module_size = sum(in_module),
      Smoking_DEG = a,
      Non_DEG = b,
      DEG_fraction = a / sum(in_module),
      OddsRatio = unname(ft$estimate),
      P = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
)

module_enrichment <- do.call(
  rbind,
  enrichment_list
)

module_enrichment$FDR <- p.adjust(
  module_enrichment$P,
  method = "BH"
)

module_enrichment <- module_enrichment[
  order(
    module_enrichment$FDR,
    module_enrichment$P
  ),
]

strong_list <- lapply(
  modules,
  function(mod) {

    in_module <- (
      gene_module$module_numeric == mod
    )

    is_deg <- gene_module$Smoking_Strong_DEG

    a <- sum(in_module & is_deg)
    b <- sum(in_module & !is_deg)

    c <- sum(!in_module & is_deg)
    d <- sum(!in_module & !is_deg)

    fisher_mat <- matrix(
      c(a, b, c, d),
      nrow = 2,
      byrow = TRUE
    )

    ft <- fisher.test(
      fisher_mat,
      alternative = "greater"
    )

    data.frame(
      module_numeric = mod,
      Strong_DEG = a,
      Strong_DEG_fraction = a / sum(in_module),
      Strong_OddsRatio = unname(ft$estimate),
      Strong_P = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
)

strong_enrichment <- do.call(
  rbind,
  strong_list
)

strong_enrichment$Strong_FDR <- p.adjust(
  strong_enrichment$Strong_P,
  method = "BH"
)

module_enrichment_final <- merge(
  module_enrichment,
  strong_enrichment,
  by = "module_numeric",
  all.x = TRUE
)

module_enrichment_final <- module_enrichment_final[
  order(
    module_enrichment_final$FDR,
    module_enrichment_final$P
  ),
]

saveRDS(
  gene_module,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_Gene_DEG_Integration.rds"
  )
)

write.table(
  gene_module,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_Gene_DEG_Integration.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  module_enrichment_final,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_Module_DEG_Enrichment.rds"
  )
)

write.table(
  module_enrichment_final,
  file.path(
    proj_dir,
    "Part2_Smoking_WGCNA_Module_DEG_Enrichment.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 9C finished\n")
cat("========================================\n")

print(module_enrichment_final)
