rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

set.seed(20260815)

required_pkgs <- c(
  "limma",
  "AnnotationDbi",
  "hugene20sttranscriptcluster.db",
  "org.Hs.eg.db",
  "ggplot2",
  "pheatmap"
)

pkg_check <- sapply(
  required_pkgs,
  requireNamespace,
  quietly = TRUE
)

print(pkg_check)

if (!all(pkg_check)) {

  stop(
    paste0(
      "Missing package(s): ",
      paste(
        names(pkg_check)[!pkg_check],
        collapse = ", "
      )
    )
  )
}

suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(hugene20sttranscriptcluster.db)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pheatmap)
})

project_dir <- file.path(project_root, "Part3_Cisplatin")

file197 <- file.path(
  project_dir,
  "GSE197561",
  "GSE197561_series_matrix.txt.gz"
)

file111 <- file.path(
  project_dir,
  "GSE111585",
  "GSE111585_series_matrix.txt.gz"
)

stopifnot(
  file.exists(file197),
  file.exists(file111)
)

out_dir <- file.path(
  project_dir,
  "Step1_FINAL_V4_OFFLINE"
)

out197 <- file.path(
  out_dir,
  "GSE197561"
)

out111 <- file.path(
  out_dir,
  "GSE111585"
)

dir.create(
  out197,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  out111,
  recursive = TRUE,
  showWarnings = FALSE
)

read_series_matrix_local <- function(file) {

  cat(
    "\nReading local file:\n",
    file,
    "\n"
  )

  con <- gzfile(
    file,
    open = "rt"
  )

  txt <- readLines(
    con,
    warn = FALSE
  )

  close(con)

  table_begin <- which(
    trimws(txt) ==
      "!series_matrix_table_begin"
  )

  table_end <- which(
    trimws(txt) ==
      "!series_matrix_table_end"
  )

  if (
    length(table_begin) != 1 ||
    length(table_end) != 1
  ) {

    stop(
      "Could not locate Series Matrix expression table."
    )
  }

  matrix_lines <- txt[
    (table_begin + 1):
      (table_end - 1)
  ]

  tc <- textConnection(
    paste(
      matrix_lines,
      collapse = "\n"
    )
  )

  dat <- read.delim(
    tc,
    header = TRUE,
    check.names = FALSE,
    quote = "\"",
    stringsAsFactors = FALSE
  )

  close(tc)

  feature_id <- as.character(
    dat[[1]]
  )

  expr <- as.matrix(
    dat[
      ,
      -1,
      drop = FALSE
    ]
  )

  storage.mode(expr) <- "numeric"

  rownames(expr) <- feature_id

  get_meta <- function(key) {

    idx <- which(
      startsWith(
        txt,
        key
      )
    )

    if (length(idx) == 0) {
      return(NULL)
    }

    one_line <- txt[
      idx[1]
    ]

    parts <- strsplit(
      one_line,
      "\t",
      fixed = TRUE
    )[[1]]

    if (length(parts) <= 1) {
      return(NULL)
    }

    values <- parts[-1]

    values <- gsub(
      '^"|"$',
      "",
      values
    )

    return(values)
  }

  gsm <- get_meta(
    "!Sample_geo_accession"
  )

  title <- get_meta(
    "!Sample_title"
  )

  platform <- get_meta(
    "!Sample_platform_id"
  )

  source <- get_meta(
    "!Sample_source_name_ch1"
  )

  data_processing <- get_meta(
    "!Sample_data_processing"
  )

  if (!is.null(gsm)) {

    if (
      length(gsm) != ncol(expr)
    ) {

      stop(
        "GSM count does not match expression columns."
      )
    }

    if (
      !all(
        colnames(expr) == gsm
      )
    ) {

      stop(
        "Expression column names do not match GSM metadata."
      )
    }
  }

  metadata <- data.frame(

    GSM = gsm,

    Title = if (!is.null(title)) {
      title
    } else {
      rep(
        NA_character_,
        ncol(expr)
      )
    },

    Source = if (!is.null(source)) {
      source
    } else {
      rep(
        NA_character_,
        ncol(expr)
      )
    },

    Platform = if (!is.null(platform)) {
      platform
    } else {
      rep(
        NA_character_,
        ncol(expr)
      )
    },

    stringsAsFactors = FALSE
  )

  return(
    list(
      expr = expr,
      metadata = metadata,
      data_processing = data_processing
    )
  )
}

clean_symbol <- function(x) {

  x <- as.character(x)

  x <- trimws(x)

  x[
    x %in% c(
      "",
      "---",
      "NA",
      "N/A",
      "NULL",
      "null"
    )
  ] <- NA

  x <- sub(
    "\\s*///.*$",
    "",
    x
  )

  x <- sub(
    "\\s*//.*$",
    "",
    x
  )

  x <- sub(
    "\\s*;.*$",
    "",
    x
  )

  x <- sub(
    "\\s*\\|.*$",
    "",
    x
  )

  x <- trimws(x)

  x[
    x == ""
  ] <- NA

  return(x)
}

collapse_gene_mean <- function(
    expr,
    symbols
) {

  keep <- !is.na(symbols) &
    symbols != ""

  expr2 <- expr[
    keep,
    ,
    drop = FALSE
  ]

  symbol2 <- symbols[
    keep
  ]

  expr_gene <- limma::avereps(
    expr2,
    ID = symbol2
  )

  good <- apply(
    expr_gene,
    1,
    function(z) {

      all(
        is.finite(z)
      ) &&
        sd(
          z,
          na.rm = TRUE
        ) > 0
    }
  )

  expr_gene <- expr_gene[
    good,
    ,
    drop = FALSE
  ]

  return(expr_gene)
}

run_qc_deg <- function(
    expr_gene,
    group,
    dataset,
    output_dir
) {

  group <- factor(
    group,
    levels = c(
      "Parental",
      "Resistant"
    )
  )

  names(group) <- colnames(
    expr_gene
  )

  cor_mat <- cor(
    expr_gene,
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  anno_col <- data.frame(
    Group = group
  )

  rownames(anno_col) <-
    names(group)

  pdf(
    file.path(
      output_dir,
      paste0(
        dataset,
        "_SampleCorrelation.pdf"
      )
    ),
    width = 7,
    height = 6
  )

  pheatmap::pheatmap(
    cor_mat,
    annotation_col = anno_col,
    annotation_row = anno_col,
    border_color = NA,
    main = paste0(
      dataset,
      " sample correlation"
    )
  )

  dev.off()

  gene_var <- apply(
    expr_gene,
    1,
    var,
    na.rm = TRUE
  )

  gene_var <- gene_var[
    is.finite(gene_var)
  ]

  top_n <- min(
    2000,
    length(gene_var)
  )

  top_genes <- names(
    sort(
      gene_var,
      decreasing = TRUE
    )
  )[seq_len(top_n)]

  pca <- prcomp(
    t(
      expr_gene[
        top_genes,
        ,
        drop = FALSE
      ]
    ),
    center = TRUE,
    scale. = TRUE
  )

  pca_percent <- (
    pca$sdev^2 /
      sum(
        pca$sdev^2
      )
  ) * 100

  pca_df <- data.frame(

    Sample = rownames(
      pca$x
    ),

    PC1 = pca$x[, 1],

    PC2 = pca$x[, 2],

    Group = group[
      rownames(
        pca$x
      )
    ],

    stringsAsFactors = FALSE
  )

  write.table(
    pca_df,
    file = file.path(
      output_dir,
      paste0(
        dataset,
        "_PCA_coordinates.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p_pca <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = Group,
      label = Sample
    )
  ) +
    geom_point(
      size = 4
    ) +
    geom_text(
      vjust = -0.8,
      size = 3
    ) +
    theme_classic(
      base_size = 12
    ) +
    labs(

      title = paste0(
        dataset,
        " PCA"
      ),

      x = sprintf(
        "PC1 (%.1f%%)",
        pca_percent[1]
      ),

      y = sprintf(
        "PC2 (%.1f%%)",
        pca_percent[2]
      )
    )

  ggsave(
    file.path(
      output_dir,
      paste0(
        dataset,
        "_PCA.png"
      )
    ),
    p_pca,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )

  ggsave(
    file.path(
      output_dir,
      paste0(
        dataset,
        "_PCA.pdf"
      )
    ),
    p_pca,
    width = 6.5,
    height = 5.5
  )

  design <- model.matrix(
    ~ 0 + group
  )

  colnames(design) <- c(
    "Parental",
    "Resistant"
  )

  contrast <- limma::makeContrasts(

    Resistant_vs_Parental =
      Resistant - Parental,

    levels = design
  )

  fit <- limma::lmFit(
    expr_gene,
    design
  )

  fit <- limma::contrasts.fit(
    fit,
    contrast
  )

  fit <- limma::eBayes(
    fit
  )

  deg <- limma::topTable(

    fit,

    coef = 1,

    number = Inf,

    adjust.method = "BH",

    sort.by = "P"
  )

  deg$Gene <- rownames(
    deg
  )

  deg <- deg[
    ,
    c(
      "Gene",
      setdiff(
        colnames(deg),
        "Gene"
      )
    )
  ]

  deg$Direction <- "NS"

  deg$Direction[
    deg$adj.P.Val < 0.05 &
      deg$logFC > 1
  ] <- "Up_in_Resistant"

  deg$Direction[
    deg$adj.P.Val < 0.05 &
      deg$logFC < -1
  ] <- "Down_in_Resistant"

  deg_sig <- deg[
    deg$Direction != "NS",
    ,
    drop = FALSE
  ]

  write.table(
    deg,
    file = file.path(
      output_dir,
      paste0(
        dataset,
        "_limma_ALL.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    deg_sig,
    file = file.path(
      output_dir,
      paste0(
        dataset,
        "_DEG_FDR005_absLogFC1.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  volcano_df <- deg

  volcano_df$minusLog10FDR <-
    -log10(
      pmax(
        volcano_df$adj.P.Val,
        1e-300
      )
    )

  p_volcano <- ggplot(
    volcano_df,
    aes(
      x = logFC,
      y = minusLog10FDR,
      color = Direction
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 1.2
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = 2
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = 2
    ) +
    theme_classic(
      base_size = 12
    ) +
    labs(
      title = paste0(
        dataset,
        ": Resistant vs Parental"
      ),
      x = "log2FC (Resistant - Parental)",
      y = "-log10(FDR)"
    )

  ggsave(
    file.path(
      output_dir,
      paste0(
        dataset,
        "_Volcano.png"
      )
    ),
    p_volcano,
    width = 6.5,
    height = 5.5,
    dpi = 300
  )

  ggsave(
    file.path(
      output_dir,
      paste0(
        dataset,
        "_Volcano.pdf"
      )
    ),
    p_volcano,
    width = 6.5,
    height = 5.5
  )

  if (
    nrow(deg_sig) >= 10
  ) {

    top50 <- head(
      deg_sig$Gene[
        order(
          deg_sig$adj.P.Val,
          -abs(
            deg_sig$logFC
          )
        )
      ],
      50
    )

  } else {

    top50 <- head(
      deg$Gene[
        order(
          deg$P.Value
        )
      ],
      50
    )
  }

  heat_mat <- expr_gene[
    top50,
    ,
    drop = FALSE
  ]

  heat_z <- t(
    scale(
      t(
        heat_mat
      )
    )
  )

  heat_z[
    !is.finite(
      heat_z
    )
  ] <- 0

  pdf(
    file.path(
      output_dir,
      paste0(
        dataset,
        "_Top50_Heatmap.pdf"
      )
    ),
    width = 7,
    height = 9
  )

  pheatmap::pheatmap(
    heat_z,
    annotation_col = anno_col,
    border_color = NA,
    main = paste0(
      dataset,
      " Top 50 genes"
    )
  )

  dev.off()

  within_cor <- c()
  between_cor <- c()

  g <- as.character(
    group
  )

  for (
    i in 1:(length(g) - 1)
  ) {

    for (
      j in (i + 1):length(g)
    ) {

      if (
        g[i] == g[j]
      ) {

        within_cor <- c(
          within_cor,
          cor_mat[i, j]
        )

      } else {

        between_cor <- c(
          between_cor,
          cor_mat[i, j]
        )
      }
    }
  }

  return(
    list(

      group = group,

      correlation = cor_mat,

      mean_within_cor = mean(
        within_cor,
        na.rm = TRUE
      ),

      mean_between_cor = mean(
        between_cor,
        na.rm = TRUE
      ),

      pca = pca_df,

      pca_percent = pca_percent,

      deg_all = deg,

      deg_sig = deg_sig,

      fit = fit
    )
  )
}

cat("\n")
cat("====================================================\n")
cat("PROCESSING GSE197561\n")
cat("====================================================\n")

geo197 <- read_series_matrix_local(
  file197
)

expr197 <- geo197$expr

meta197 <- geo197$metadata

cat(
  "Matrix:",
  nrow(expr197),
  "features x",
  ncol(expr197),
  "samples\n"
)

print(
  meta197
)

group_map197 <- c(

  GSM5920580 = "Resistant",
  GSM5920581 = "Resistant",
  GSM5920582 = "Resistant",

  GSM5920583 = "Parental",
  GSM5920584 = "Parental",
  GSM5920585 = "Parental"
)

gsm197 <- colnames(
  expr197
)

stopifnot(
  all(
    gsm197 %in%
      names(group_map197)
  )
)

group197 <- factor(

  group_map197[
    gsm197
  ],

  levels = c(
    "Parental",
    "Resistant"
  )
)

names(group197) <- gsm197

audit197 <- meta197

audit197$Group <-
  as.character(
    group197
  )

cat(
  "\nGSE197561 sample audit:\n"
)

print(
  audit197
)

stopifnot(
  sum(
    group197 ==
      "Parental"
  ) == 3,

  sum(
    group197 ==
      "Resistant"
  ) == 3
)

q197 <- quantile(

  as.numeric(
    expr197
  ),

  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.99,
    1
  ),

  na.rm = TRUE
)

cat(
  "\nGSE197561 quantiles:\n"
)

print(
  q197
)

if (
  q197[5] > 100
) {

  stop(
    "GSE197561 unexpectedly looks non-log2."
  )
}

cat(
  "GSE197561: already log2/RMA scale. No transformation.\n"
)

probe197 <- as.character(
  rownames(
    expr197
  )
)

cat(
  "\nMapping GSE197561 feature IDs...\n"
)

map197 <- suppressMessages(

  AnnotationDbi::select(

    hugene20sttranscriptcluster.db,

    keys = probe197,

    columns = "SYMBOL",

    keytype = "PROBEID"
  )
)

map197 <- map197[
  !is.na(
    map197$SYMBOL
  ),
  ,
  drop = FALSE
]

map197 <- map197[
  !duplicated(
    map197$PROBEID
  ),
  ,
  drop = FALSE
]

symbol197 <- map197$SYMBOL[
  match(
    probe197,
    map197$PROBEID
  )
]

symbol197 <- clean_symbol(
  symbol197
)

mapped197 <- mean(
  !is.na(
    symbol197
  )
) * 100

cat(
  "Mapped percentage:",
  round(
    mapped197,
    2
  ),
  "%\n"
)

if (
  mapped197 < 20
) {

  stop(
    "GSE197561 mapping unexpectedly low."
  )
}

expr197_gene <- collapse_gene_mean(

  expr = expr197,

  symbols = symbol197
)

cat(
  "GSE197561 gene-level matrix:",
  nrow(expr197_gene),
  "genes x",
  ncol(expr197_gene),
  "samples\n"
)

saveRDS(

  expr197_gene,

  file.path(
    out197,
    "GSE197561_GeneLevel_Expression.rds"
  )
)

write.table(

  data.frame(
    Feature_ID = probe197,
    GeneSymbol = symbol197
  ),

  file = file.path(
    out197,
    "GSE197561_Feature_to_Gene.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)

res197 <- run_qc_deg(

  expr_gene = expr197_gene,

  group = group197,

  dataset = "GSE197561",

  output_dir = out197
)

cat("\n")
cat("====================================================\n")
cat("PROCESSING GSE111585\n")
cat("====================================================\n")

geo111 <- read_series_matrix_local(
  file111
)

expr111_raw <- geo111$expr

meta111 <- geo111$metadata

cat(
  "Matrix:",
  nrow(expr111_raw),
  "features x",
  ncol(expr111_raw),
  "samples\n"
)

print(
  meta111
)

group_map111 <- c(

  GSM3034616 = "Parental",
  GSM3034617 = "Parental",
  GSM3034618 = "Parental",

  GSM3034619 = "Resistant",
  GSM3034620 = "Resistant",
  GSM3034621 = "Resistant"
)

gsm111 <- colnames(
  expr111_raw
)

stopifnot(
  all(
    gsm111 %in%
      names(group_map111)
  )
)

group111 <- factor(

  group_map111[
    gsm111
  ],

  levels = c(
    "Parental",
    "Resistant"
  )
)

names(group111) <- gsm111

audit111 <- meta111

audit111$Group <-
  as.character(
    group111
  )

cat(
  "\nGSE111585 sample audit:\n"
)

print(
  audit111
)

stopifnot(
  sum(
    group111 ==
      "Parental"
  ) == 3,

  sum(
    group111 ==
      "Resistant"
  ) == 3
)

q111_before <- quantile(

  as.numeric(
    expr111_raw
  ),

  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.99,
    1
  ),

  na.rm = TRUE
)

cat(
  "\nGSE111585 before log2:\n"
)

print(
  q111_before
)

expr111 <- log2(
  pmax(
    expr111_raw,
    1
  )
)

q111_after <- quantile(

  as.numeric(
    expr111
  ),

  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.99,
    1
  ),

  na.rm = TRUE
)

cat(
  "\nGSE111585 after log2:\n"
)

print(
  q111_after
)

probe111 <- as.character(
  rownames(
    expr111
  )
)

cat(
  "\nMapping GSE111585 GenBank accessions...\n"
)

map111 <- suppressMessages(

  AnnotationDbi::select(

    org.Hs.eg.db,

    keys = probe111,

    columns = "SYMBOL",

    keytype = "ACCNUM"
  )
)

map111 <- map111[
  !is.na(
    map111$SYMBOL
  ),
  ,
  drop = FALSE
]

map111 <- map111[
  !duplicated(
    map111$ACCNUM
  ),
  ,
  drop = FALSE
]

symbol111 <- map111$SYMBOL[
  match(
    probe111,
    map111$ACCNUM
  )
]

symbol111 <- clean_symbol(
  symbol111
)

mapped111 <- mean(
  !is.na(
    symbol111
  )
) * 100

cat(
  "Mapped percentage:",
  round(
    mapped111,
    2
  ),
  "%\n"
)

cat(
  "Mapped features:",
  sum(
    !is.na(
      symbol111
    )
  ),
  "/",
  length(
    symbol111
  ),
  "\n"
)

if (
  sum(
    !is.na(
      symbol111
    )
  ) < 5000
) {

  stop(
    "Too few GSE111585 coding genes mapped. Stop for annotation review."
  )
}

expr111_gene <- collapse_gene_mean(

  expr = expr111,

  symbols = symbol111
)

cat(
  "GSE111585 gene-level matrix:",
  nrow(expr111_gene),
  "genes x",
  ncol(expr111_gene),
  "samples\n"
)

saveRDS(

  expr111_gene,

  file.path(
    out111,
    "GSE111585_GeneLevel_Expression.rds"
  )
)

write.table(

  data.frame(
    Feature_ID = probe111,
    GeneSymbol = symbol111
  ),

  file = file.path(
    out111,
    "GSE111585_Feature_to_Gene.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)

res111 <- run_qc_deg(

  expr_gene = expr111_gene,

  group = group111,

  dataset = "GSE111585",

  output_dir = out111
)

summary_all <- data.frame(

  Dataset = c(
    "GSE197561",
    "GSE111585"
  ),

  Role = c(
    "Discovery",
    "Independent_replication"
  ),

  Cell_line = c(
    "UMSCC1",
    "SCC9"
  ),

  N_Parental = c(
    3,
    3
  ),

  N_Resistant = c(
    3,
    3
  ),

  Original_features = c(
    nrow(
      expr197
    ),
    nrow(
      expr111_raw
    )
  ),

  Gene_level_features = c(
    nrow(
      expr197_gene
    ),
    nrow(
      expr111_gene
    )
  ),

  Mapped_percent = c(
    round(
      mapped197,
      2
    ),
    round(
      mapped111,
      2
    )
  ),

  PC1_percent = c(
    round(
      res197$pca_percent[1],
      2
    ),
    round(
      res111$pca_percent[1],
      2
    )
  ),

  PC2_percent = c(
    round(
      res197$pca_percent[2],
      2
    ),
    round(
      res111$pca_percent[2],
      2
    )
  ),

  Mean_within_group_cor = c(
    round(
      res197$mean_within_cor,
      4
    ),
    round(
      res111$mean_within_cor,
      4
    )
  ),

  Mean_between_group_cor = c(
    round(
      res197$mean_between_cor,
      4
    ),
    round(
      res111$mean_between_cor,
      4
    )
  ),

  DEG_total = c(
    nrow(
      res197$deg_sig
    ),
    nrow(
      res111$deg_sig
    )
  ),

  DEG_up = c(

    sum(
      res197$deg_all$Direction ==
        "Up_in_Resistant"
    ),

    sum(
      res111$deg_all$Direction ==
        "Up_in_Resistant"
    )
  ),

  DEG_down = c(

    sum(
      res197$deg_all$Direction ==
        "Down_in_Resistant"
    ),

    sum(
      res111$deg_all$Direction ==
        "Down_in_Resistant"
    )
  ),

  stringsAsFactors = FALSE
)

write.table(

  summary_all,

  file = file.path(
    out_dir,
    "Part3_Step1_FINAL_V4_Summary.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)

result197 <- list(

  Dataset = "GSE197561",

  Role = "Discovery",

  audit = audit197,

  group = group197,

  expression_gene = expr197_gene,

  mapped_percent = mapped197,

  pca = res197$pca,

  pca_percent = res197$pca_percent,

  correlation = res197$correlation,

  deg_all = res197$deg_all,

  deg_sig = res197$deg_sig
)

result111 <- list(

  Dataset = "GSE111585",

  Role = "Independent_replication",

  audit = audit111,

  group = group111,

  expression_gene = expr111_gene,

  mapped_percent = mapped111,

  pca = res111$pca,

  pca_percent = res111$pca_percent,

  correlation = res111$correlation,

  deg_all = res111$deg_all,

  deg_sig = res111$deg_sig
)

results <- list(

  GSE197561 = result197,

  GSE111585 = result111
)

saveRDS(

  results,

  file = file.path(
    out_dir,
    "Part3_Step1_FINAL_V4_AllResults.rds"
  )
)

capture.output(

  sessionInfo(),

  file = file.path(
    out_dir,
    "Part3_Step1_FINAL_V4_sessionInfo.txt"
  )
)

cat("\n\n")
cat("####################################################\n")
cat("STEP 1 FINAL V4 FINISHED SUCCESSFULLY\n")
cat("####################################################\n\n")

print(
  summary_all
)

cat(
  "\nResults saved to:\n",
  out_dir,
  "\n"
)
