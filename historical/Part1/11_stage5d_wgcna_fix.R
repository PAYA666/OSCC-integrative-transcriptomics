rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(WGCNA)
library(survival)

options(stringsAsFactors = FALSE)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

moduleColors <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_moduleColors.rds"
)

tumor_meta <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_tumor_metadata.rds"
)

clinical_model <- readRDS(
  "TCGA_OSCC_Nature8_clinical_model_raw.rds"
)

ME_obj <- moduleEigengenes(
  datExpr,
  colors = moduleColors
)

MEs_color <- orderMEs(
  ME_obj$eigengenes
)

cat("Module eigengenes:\n")
print(colnames(MEs_color))

MEs_use <- MEs_color[
  ,
  colnames(MEs_color) != "MEgrey",
  drop = FALSE
]

cat("\nBiological modules:\n")
print(colnames(MEs_use))

sample_map <- tumor_meta[
  match(
    rownames(datExpr),
    tumor_meta$sample_barcode
  ),
  c("sample_barcode", "patient_barcode")
]

clinical_wgcna <- clinical_model[
  match(
    sample_map$patient_barcode,
    clinical_model$patient_barcode
  ),
]

stopifnot(
  identical(
    clinical_wgcna$patient_barcode,
    sample_map$patient_barcode
  )
)

sex_num <- ifelse(
  clinical_wgcna$sex == "male", 1,
  ifelse(clinical_wgcna$sex == "female", 0, NA)
)

stage_num <- c(
  "Stage I" = 1,
  "Stage II" = 2,
  "Stage III" = 3,
  "Stage IV" = 4,
  "Stage IVA" = 4,
  "Stage IVB" = 4,
  "Stage IVC" = 4
)[clinical_wgcna$stage]

T_num <- ifelse(
  clinical_wgcna$T == "T1", 1,
  ifelse(
    clinical_wgcna$T == "T2", 2,
    ifelse(
      clinical_wgcna$T == "T3", 3,
      ifelse(
        clinical_wgcna$T %in% c("T4","T4a","T4b"), 4, NA
      )
    )
  )
)

N_num <- ifelse(
  clinical_wgcna$N == "N0", 0,
  ifelse(
    clinical_wgcna$N == "N1", 1,
    ifelse(
      clinical_wgcna$N %in% c("N2","N2a","N2b","N2c"), 2,
      ifelse(clinical_wgcna$N == "N3", 3, NA)
    )
  )
)

traitData <- data.frame(
  Age   = clinical_wgcna$age,
  Sex   = sex_num,
  Stage = as.numeric(stage_num),
  T     = T_num,
  N     = N_num,
  Death = clinical_wgcna$OS_event
)

rownames(traitData) <- rownames(datExpr)

moduleTraitCor <- cor(
  MEs_use,
  traitData,
  use = "pairwise.complete.obs"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(datExpr)
)

textMatrix <- paste0(
  sprintf("%.2f", moduleTraitCor),
  "\n(",
  sprintf("%.2g", moduleTraitPvalue),
  ")"
)

dim(textMatrix) <- dim(moduleTraitCor)

pdf(
  "Figure_Nature8_OSCC_WGCNA_ModuleTraitHeatmap_COLOR.pdf",
  width = 10,
  height = 7
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),
  yLabels = colnames(MEs_use),
  ySymbols = colnames(MEs_use),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1,1),
  main = "WGCNA module-trait relationships"
)

dev.off()

surv_df <- data.frame(
  OS_time = clinical_wgcna$OS_time,
  OS_event = clinical_wgcna$OS_event,
  MEs_use,
  check.names = FALSE
)

surv_df <- surv_df[
  !is.na(surv_df$OS_time) &
    surv_df$OS_time > 0 &
    !is.na(surv_df$OS_event),
]

cox_list <- lapply(
  colnames(MEs_use),
  function(me) {

    fit <- coxph(
      as.formula(
        paste0(
          "Surv(OS_time, OS_event) ~ `",
          me,
          "`"
        )
      ),
      data = surv_df
    )

    s <- summary(fit)

    data.frame(
      Module = me,
      HR = s$coefficients[1,"exp(coef)"],
      lower95 = s$conf.int[1,"lower .95"],
      upper95 = s$conf.int[1,"upper .95"],
      P = s$coefficients[1,"Pr(>|z|)"]
    )
  }
)

module_cox_color <- do.call(
  rbind,
  cox_list
)

module_cox_color$FDR <- p.adjust(
  module_cox_color$P,
  method = "BH"
)

module_cox_color <- module_cox_color[
  order(module_cox_color$P),
]

saveRDS(
  MEs_color,
  "TCGA_OSCC_Nature8_WGCNA_moduleEigengenes_COLOR.rds"
)

saveRDS(
  module_cox_color,
  "TCGA_OSCC_Nature8_WGCNA_Module_OS_Cox_COLOR.rds"
)

write.table(
  module_cox_color,
  "TCGA_OSCC_Nature8_WGCNA_Module_OS_Cox_COLOR.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  moduleTraitCor,
  "TCGA_OSCC_Nature8_WGCNA_ModuleTrait_Correlation_COLOR.tsv",
  sep = "\t",
  quote = FALSE
)

write.table(
  moduleTraitPvalue,
  "TCGA_OSCC_Nature8_WGCNA_ModuleTrait_Pvalue_COLOR.tsv",
  sep = "\t",
  quote = FALSE
)
