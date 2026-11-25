project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

proj_dir <- project_root

expr_long <- readRDS(
  file.path(
    proj_dir,
    "Part2_Smoking_KeyGenes_VST_LongFormat.rds"
  )
)

exposure <- readRDS(
  file.path(
    proj_dir,
    "Part2_TCGA_OSCC_Nature8_Exposure_341_CLEAN.rds"
  )
)

stopifnot(
  length(unique(exposure$patient_id)) ==
    nrow(exposure)
)

cat("\nOriginal smoking categories:\n")

print(
  table(
    exposure$tobacco_smoking_status,
    useNA = "ifany"
  )
)

status_raw <- trimws(
  as.character(
    exposure$tobacco_smoking_status
  )
)

status_lower <- tolower(
  status_raw
)

exposure$smoking3 <- NA_character_

exposure$smoking3[
  !is.na(status_lower) &
  (
    grepl("lifelong", status_lower) |
    grepl("never", status_lower)
  )
] <- "Never"

exposure$smoking3[
  !is.na(status_lower) &
  grepl("current", status_lower)
] <- "Current"

exposure$smoking3[
  !is.na(status_lower) &
  (
    grepl("reformed", status_lower) |
    grepl("former", status_lower)
  )
] <- "Former"

unmapped <- unique(
  status_raw[
    !is.na(status_raw) &
    status_raw != "" &
    is.na(exposure$smoking3)
  ]
)

if (length(unmapped) > 0) {

  cat("\nUnmapped Smoking categories:\n")
  print(unmapped)

  stop(
    "There are unmapped Smoking categories. Check before continuing."
  )
}

exposure$smoking3 <- factor(
  exposure$smoking3,
  levels = c(
    "Never",
    "Former",
    "Current"
  )
)

cat("\nCollapsed Smoking groups:\n")

print(
  table(
    exposure$smoking3,
    useNA = "ifany"
  )
)

exposure$pack_years_numeric <-
  suppressWarnings(
    as.numeric(
      as.character(
        exposure$pack_years_smoked
      )
    )
  )

cat("\nPack-years availability:\n")

cat(
  "Non-missing:",
  sum(
    !is.na(exposure$pack_years_numeric)
  ),
  "/",
  nrow(exposure),
  "\n"
)

cat(
  "Positive pack-years:",
  sum(
    exposure$pack_years_numeric > 0,
    na.rm = TRUE
  ),
  "\n"
)

expo_sub <- exposure[
  ,
  c(
    "patient_id",
    "smoking3",
    "smoking_bin",
    "pack_years_numeric"
  )
]

dat <- merge(
  expr_long,
  expo_sub,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c(
    "_expr",
    "_exposure"
  )
)

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

stopifnot(
  all(
    key_genes %in%
      unique(dat$gene_name)
  )
)

status_summary <- do.call(
  rbind,
  lapply(
    key_genes,
    function(g) {

      tmp <- dat[
        dat$gene_name == g &
        !is.na(dat$smoking3),
      ]

      do.call(
        rbind,
        lapply(
          levels(exposure$smoking3),
          function(grp) {

            x <- tmp$VST_expression[
              tmp$smoking3 == grp
            ]

            data.frame(
              gene_name = g,
              Smoking_status = grp,
              N = length(x),
              Mean = mean(
                x,
                na.rm = TRUE
              ),
              Median = median(
                x,
                na.rm = TRUE
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

kw_results <- do.call(
  rbind,
  lapply(
    key_genes,
    function(g) {

      tmp <- dat[
        dat$gene_name == g &
        !is.na(dat$smoking3),
      ]

      kt <- kruskal.test(
        VST_expression ~ smoking3,
        data = tmp
      )

      data.frame(
        gene_name = g,
        P = kt$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

kw_results$FDR <- p.adjust(
  kw_results$P,
  method = "BH"
)

kw_results <- kw_results[
  order(kw_results$FDR),
]

comparisons <- list(
  c("Current", "Never"),
  c("Former", "Never"),
  c("Current", "Former")
)

pairwise_results <- do.call(
  rbind,
  lapply(
    key_genes,
    function(g) {

      tmp <- dat[
        dat$gene_name == g &
        !is.na(dat$smoking3),
      ]

      do.call(
        rbind,
        lapply(
          comparisons,
          function(comp) {

            group1 <- comp[1]
            group2 <- comp[2]

            x1 <- tmp$VST_expression[
              tmp$smoking3 == group1
            ]

            x2 <- tmp$VST_expression[
              tmp$smoking3 == group2
            ]

            wt <- wilcox.test(
              x1,
              x2,
              exact = FALSE
            )

            data.frame(
              gene_name = g,

              Comparison =
                paste0(
                  group1,
                  "_vs_",
                  group2
                ),

              Group1_Median =
                median(
                  x1,
                  na.rm = TRUE
                ),

              Group2_Median =
                median(
                  x2,
                  na.rm = TRUE
                ),

              Median_Difference =
                median(
                  x1,
                  na.rm = TRUE
                ) -
                median(
                  x2,
                  na.rm = TRUE
                ),

              P = wt$p.value,

              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

pairwise_results$FDR <- p.adjust(
  pairwise_results$P,
  method = "BH"
)

pairwise_results <- pairwise_results[
  order(
    pairwise_results$FDR,
    pairwise_results$P
  ),
]

pack_dat <- dat[
  dat$smoking_bin_exposure == "Ever" &
  !is.na(dat$pack_years_numeric) &
  dat$pack_years_numeric > 0,
]

packyear_results <- do.call(
  rbind,
  lapply(
    key_genes,
    function(g) {

      tmp <- pack_dat[
        pack_dat$gene_name == g,
      ]

      ct <- cor.test(
        tmp$VST_expression,
        tmp$pack_years_numeric,
        method = "spearman",
        exact = FALSE
      )

      data.frame(
        gene_name = g,

        N =
          nrow(tmp),

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
)

packyear_results$FDR <- p.adjust(
  packyear_results$P,
  method = "BH"
)

packyear_results <- packyear_results[
  order(
    packyear_results$FDR,
    packyear_results$P
  ),
]

write.table(
  status_summary,
  file.path(
    proj_dir,
    "Part2_Smoking3_KeyGenes_ExpressionSummary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  kw_results,
  file.path(
    proj_dir,
    "Part2_Smoking3_KeyGenes_KruskalWallis.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pairwise_results,
  file.path(
    proj_dir,
    "Part2_Smoking3_KeyGenes_Pairwise.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  packyear_results,
  file.path(
    proj_dir,
    "Part2_PackYears_KeyGenes_Spearman.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  list(
    status_summary = status_summary,
    kw_results = kw_results,
    pairwise_results = pairwise_results,
    packyear_results = packyear_results
  ),
  file.path(
    proj_dir,
    "Part2_Smoking_ExposureDose_Sensitivity_FINAL.rds"
  )
)

cat("\n=====================================\n")
cat("Step 14 finished\n")
cat("=====================================\n")

cat("\nSmoking3:\n")
print(
  table(
    exposure$smoking3,
    useNA = "ifany"
  )
)

cat("\nGlobal Smoking-status tests:\n")
print(kw_results)

cat("\nPack-years correlations:\n")
print(packyear_results)
