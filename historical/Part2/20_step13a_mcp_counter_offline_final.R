project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

stopifnot(
  exists("expr_mcp"),
  ncol(expr_mcp) == 333
)

cat(
  "Expression matrix:",
  nrow(expr_mcp),
  "genes x",
  ncol(expr_mcp),
  "samples\n"
)

mcp_markers <- list(

  "T cells" = c(
    "CD28", "CD3D", "CD3G", "CD5", "CD6",
    "CHRM3-AS2", "CTLA4", "FLT3LG", "ICOS",
    "MAL", "MGC40069", "PBX4", "SIRPG",
    "THEMIS", "TNFRSF25", "TRAT1"
  ),

  "CD8 T cells" = c(
    "CD8B"
  ),

  "Cytotoxic lymphocytes" = c(
    "CD8A", "EOMES", "FGFBP2", "GNLY",
    "KLRC3", "KLRC4", "KLRD1"
  ),

  "B lineage" = c(
    "BANK1", "CD19", "CD22", "CD79A",
    "CR2", "FCRL2", "IGKC", "MS4A1",
    "PAX5"
  ),

  "NK cells" = c(
    "CD160", "KIR2DL1", "KIR2DL3",
    "KIR2DL4", "KIR3DL1", "KIR3DS1",
    "NCR1", "PTGDR", "SH2D1B"
  ),

  "Monocytic lineage" = c(
    "ADAP2", "CSF1R", "FPR3", "KYNU",
    "PLA2G7", "RASSF4", "TFEC"
  ),

  "Myeloid dendritic cells" = c(
    "CD1A", "CD1B", "CD1E", "CLEC10A",
    "CLIC2", "WFDC21P"
  ),

  "Neutrophils" = c(
    "CA4", "CEACAM3", "CXCR1", "CXCR2",
    "CYP4F3", "FCGR3B", "HAL", "KCNJ15",
    "MEGF9", "SLC25A37", "STEAP4",
    "TECPR2", "TLE3", "TNFRSF10C", "VNN3"
  ),

  "Endothelial cells" = c(
    "ACVRL1", "APLN", "BCL6B", "BMP6",
    "BMX", "CDH5", "CLEC14A", "DIPK2B",
    "EDN1", "ADGRL4", "EMCN", "ESAM",
    "ESM1", "FAM124B", "HECW2", "HHIP",
    "KDR", "MMRN1", "MMRN2", "MYCT1",
    "PALMD", "PEAR1", "PGF", "PLXNA2",
    "PTPRB", "ROBO4", "C", "SHANK3",
    "SHE", "TEK", "TIE1", "VEPH1",
    "VWF"
  ),

  "Fibroblasts" = c(
    "COL1A1", "COL3A1", "COL6A1",
    "COL6A2", "DCN", "GREM1",
    "PAMR1", "TAGLN"
  )
)

mcp_markers_present <- lapply(
  mcp_markers,
  function(x) {
    intersect(
      x,
      rownames(expr_mcp)
    )
  }
)

marker_coverage <- data.frame(

  CellType = names(mcp_markers),

  Expected_markers =
    lengths(mcp_markers),

  Available_markers =
    lengths(mcp_markers_present),

  Coverage =
    lengths(mcp_markers_present) /
    lengths(mcp_markers),

  stringsAsFactors = FALSE
)

cat("\nMarker coverage:\n")

print(marker_coverage)

if (any(
  marker_coverage$Available_markers == 0
)) {

  stop(
    paste(
      "No markers detected for:",
      paste(
        marker_coverage$CellType[
          marker_coverage$Available_markers == 0
        ],
        collapse = ", "
      )
    )
  )
}

mcp_scores <- do.call(
  rbind,
  lapply(
    mcp_markers_present,
    function(markers) {

      colMeans(
        expr_mcp[
          markers,
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )

    }
  )
)

rownames(mcp_scores) <-
  names(mcp_markers_present)

colnames(mcp_scores) <-
  colnames(expr_mcp)

cat(
  "\nMCP scores:",
  nrow(mcp_scores),
  "populations x",
  ncol(mcp_scores),
  "samples\n"
)

mcp_df <- as.data.frame(
  t(mcp_scores)
)

mcp_df$sample_barcode <-
  rownames(mcp_df)

if (!exists("sample_info333")) {

  sample_info <- readRDS(
    file.path(
      proj_dir,
      "Part2_Smoking_sample_metadata_patient341.rds"
    )
  )

  sample_info333 <- sample_info[
    !is.na(sample_info$smoking_bin),
  ]
}

stopifnot(
  nrow(sample_info333) == 333
)

idx_meta <- match(
  mcp_df$sample_barcode,
  sample_info333$sample_barcode
)

mcp_df$patient_id <-
  sample_info333$patient_id[idx_meta]

mcp_df$smoking_bin <-
  sample_info333$smoking_bin[idx_meta]

stopifnot(
  all(!is.na(mcp_df$smoking_bin))
)

cat("\nSmoking groups:\n")

print(
  table(mcp_df$smoking_bin)
)

cell_types <- names(
  mcp_markers_present
)

tme_results <- lapply(
  cell_types,
  function(cell) {

    never_values <- as.numeric(
      mcp_df[
        mcp_df$smoking_bin == "Never",
        cell
      ]
    )

    ever_values <- as.numeric(
      mcp_df[
        mcp_df$smoking_bin == "Ever",
        cell
      ]
    )

    wt <- wilcox.test(
      ever_values,
      never_values,
      exact = FALSE
    )

    data.frame(

      CellType = cell,

      Never_Mean =
        mean(
          never_values,
          na.rm = TRUE
        ),

      Ever_Mean =
        mean(
          ever_values,
          na.rm = TRUE
        ),

      Never_Median =
        median(
          never_values,
          na.rm = TRUE
        ),

      Ever_Median =
        median(
          ever_values,
          na.rm = TRUE
        ),

      Median_Difference =
        median(
          ever_values,
          na.rm = TRUE
        ) -
        median(
          never_values,
          na.rm = TRUE
        ),

      P = wt$p.value,

      stringsAsFactors = FALSE
    )
  }
)

tme_results <- do.call(
  rbind,
  tme_results
)

tme_results$FDR <- p.adjust(
  tme_results$P,
  method = "BH"
)

tme_results$Direction <- ifelse(
  tme_results$Median_Difference > 0,
  "Ever_higher",
  "Never_higher"
)

tme_results <- tme_results[
  order(
    tme_results$FDR,
    tme_results$P
  ),
]

key_genes <- c(
  "CYP1A1",
  "CES1",
  "SFRP4",
  "SFRP2",
  "PDGFRA",
  "OGN",
  "FNDC1",
  "TK1",
  "MYBL2"
)

genes_present <- intersect(
  key_genes,
  rownames(expr_mcp)
)

fibro_score <- as.numeric(
  mcp_df[
    match(
      colnames(expr_mcp),
      mcp_df$sample_barcode
    ),
    "Fibroblasts"
  ]
)

gene_fibro_cor <- lapply(
  genes_present,
  function(g) {

    ct <- cor.test(
      as.numeric(expr_mcp[g, ]),
      fibro_score,
      method = "spearman",
      exact = FALSE
    )

    data.frame(

      gene_name = g,

      Spearman_Rho =
        unname(
          ct$estimate
        ),

      P =
        ct$p.value,

      stringsAsFactors = FALSE
    )
  }
)

gene_fibro_cor <- do.call(
  rbind,
  gene_fibro_cor
)

gene_fibro_cor$FDR <- p.adjust(
  gene_fibro_cor$P,
  method = "BH"
)

gene_fibro_cor <- gene_fibro_cor[
  order(
    -abs(
      gene_fibro_cor$Spearman_Rho
    )
  ),
]

saveRDS(
  mcp_scores,
  file.path(
    proj_dir,
    "Part2_Smoking_MCPcounter_Scores_OFFLINE.rds"
  )
)

write.table(
  marker_coverage,
  file.path(
    proj_dir,
    "Part2_Smoking_MCPcounter_MarkerCoverage.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  tme_results,
  file.path(
    proj_dir,
    "Part2_Smoking_MCPcounter_Ever_vs_Never.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_fibro_cor,
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_FibroblastCorrelation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat("Step 13A OFFLINE finished\n")
cat("========================================\n")

cat("\nTME comparison:\n")
print(tme_results)

cat("\nKey genes vs Fibroblast score:\n")
print(gene_fibro_cor)
