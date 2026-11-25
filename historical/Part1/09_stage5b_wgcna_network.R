project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(WGCNA)

setwd(project_root)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

sample_qc <- read.delim(
  "TCGA_OSCC_Nature8_WGCNA_sample_QC.tsv",
  stringsAsFactors = FALSE
)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

candidate_outliers <- sample_qc$sample[
  sample_qc$connectivity_Z < -4
]

cat("Candidate outliers (Z < -4):\n")
print(candidate_outliers)

datExpr_all <- datExpr

datExpr_trim <- datExpr[
  !rownames(datExpr) %in% candidate_outliers,
  ,
  drop = FALSE
]

cat("\nDimensions:\n")
cat("All:", dim(datExpr_all), "\n")
cat("Trimmed:", dim(datExpr_trim), "\n")

powers <- c(
  1:10,
  seq(12, 30, by = 2)
)

sft_all <- pickSoftThreshold(
  datExpr_all,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 5
)

sft_trim <- pickSoftThreshold(
  datExpr_trim,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 5
)

saveRDS(
  sft_all,
  "TCGA_OSCC_Nature8_WGCNA_softThreshold_ALL.rds"
)

saveRDS(
  sft_trim,
  "TCGA_OSCC_Nature8_WGCNA_softThreshold_TRIMMED.rds"
)

pdf(
  "Figure_Nature8_OSCC_WGCNA_SoftThreshold_Comparison.pdf",
  width = 12,
  height = 10
)

par(mfrow = c(2,2))

plot(
  sft_all$fitIndices[,1],
  -sign(sft_all$fitIndices[,3]) * sft_all$fitIndices[,2],
  xlab = "Soft Threshold (power)",
  ylab = "Signed R^2",
  main = "All samples: scale independence",
  type = "n"
)

text(
  sft_all$fitIndices[,1],
  -sign(sft_all$fitIndices[,3]) * sft_all$fitIndices[,2],
  labels = powers,
  cex = 0.8
)

abline(h = 0.80, lty = 2)
abline(h = 0.85, lty = 3)

plot(
  sft_all$fitIndices[,1],
  sft_all$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean connectivity",
  main = "All samples: mean connectivity",
  type = "n"
)

text(
  sft_all$fitIndices[,1],
  sft_all$fitIndices[,5],
  labels = powers,
  cex = 0.8
)

plot(
  sft_trim$fitIndices[,1],
  -sign(sft_trim$fitIndices[,3]) * sft_trim$fitIndices[,2],
  xlab = "Soft Threshold (power)",
  ylab = "Signed R^2",
  main = "Trimmed: scale independence",
  type = "n"
)

text(
  sft_trim$fitIndices[,1],
  -sign(sft_trim$fitIndices[,3]) * sft_trim$fitIndices[,2],
  labels = powers,
  cex = 0.8
)

abline(h = 0.80, lty = 2)
abline(h = 0.85, lty = 3)

plot(
  sft_trim$fitIndices[,1],
  sft_trim$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean connectivity",
  main = "Trimmed: mean connectivity",
  type = "n"
)

text(
  sft_trim$fitIndices[,1],
  sft_trim$fitIndices[,5],
  labels = powers,
  cex = 0.8
)

dev.off()

summary_all <- data.frame(
  power = sft_all$fitIndices[,1],
  signed_R2 =
    -sign(sft_all$fitIndices[,3]) *
    sft_all$fitIndices[,2],
  mean_connectivity =
    sft_all$fitIndices[,5]
)

summary_trim <- data.frame(
  power = sft_trim$fitIndices[,1],
  signed_R2 =
    -sign(sft_trim$fitIndices[,3]) *
    sft_trim$fitIndices[,2],
  mean_connectivity =
    sft_trim$fitIndices[,5]
)

write.table(
  summary_all,
  "TCGA_OSCC_Nature8_WGCNA_softThreshold_ALL.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_trim,
  "TCGA_OSCC_Nature8_WGCNA_softThreshold_TRIMMED.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
