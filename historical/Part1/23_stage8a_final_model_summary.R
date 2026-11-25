rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(timeROC)
library(ggplot2)
library(dplyr)

fit_z <- readRDS(
  "TCGA_OSCC_Final8Gene_STANDARDIZED_CoxModel.rds"
)

beta_z <- readRDS(
  "TCGA_OSCC_Final8Gene_STANDARDIZED_beta.rds"
)

tcga_final <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_RiskScore.rds"
)

nested_oof <- readRDS(
  "TCGA_OSCC_Nature8_Stage6F_NestedCV_PatientPredictions.rds"
)

geo_validation <- readRDS(
  "GSE41613_Final8Gene_ExternalValidation.rds"
)

combined_fit <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_CombinedClinicalModel.rds"
)

clinical_dat <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalCompleteCases.rds"
)

formula_text <- paste0(
  "Risk score = ",
  paste(
    sprintf(
      "%+.6f × Z(%s)",
      beta_z,
      names(beta_z)
    ),
    collapse = " "
  )
)

formula_text <- sub(
  "^Risk score = \\+",
  "Risk score = ",
  formula_text
)

cat("\n========================================\n")
cat("FINAL STANDARDIZED MODEL\n")
cat("========================================\n")
cat(formula_text, "\n")

writeLines(
  formula_text,
  "TCGA_OSCC_Final8Gene_STANDARDIZED_RiskScore_Formula.txt"
)

coef_formula <- data.frame(
  gene_name = names(beta_z),
  standardized_beta = as.numeric(beta_z),
  direction = ifelse(
    beta_z > 0,
    "Risk",
    "Protective"
  ),
  stringsAsFactors = FALSE
)

write.table(
  coef_formula,
  "TCGA_OSCC_Final8Gene_STANDARDIZED_Formula.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

s_z <- summary(fit_z)

forest8 <- data.frame(
  gene_name = rownames(s_z$coefficients),
  beta = s_z$coefficients[, "coef"],
  HR = s_z$coefficients[, "exp(coef)"],
  lower95 = s_z$conf.int[, "lower .95"],
  upper95 = s_z$conf.int[, "upper .95"],
  P = s_z$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

forest8$label <- paste0(
  forest8$gene_name,
  "   HR ",
  sprintf("%.2f", forest8$HR),
  " (",
  sprintf("%.2f", forest8$lower95),
  "–",
  sprintf("%.2f", forest8$upper95),
  ")"
)

forest8$gene_name <- factor(
  forest8$gene_name,
  levels = rev(forest8$gene_name)
)

p_forest <- ggplot(
  forest8,
  aes(
    x = HR,
    y = gene_name
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = lower95,
      xmax = upper95
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_log10() +
  labs(
    title = "Multivariable Cox analysis of the 8-gene signature",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 11
    )
  )

ggsave(
  "Figure_Final8Gene_Forest_PUBLICATION.pdf",
  p_forest,
  width = 8,
  height = 6
)

ggsave(
  "Figure_Final8Gene_Forest_PUBLICATION.png",
  p_forest,
  width = 8,
  height = 6,
  dpi = 600
)

write.table(
  forest8,
  "TCGA_OSCC_Final8Gene_Forest_Table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

c_tcga <- survival::concordance(
  Surv(
    OS_time,
    OS_event
  ) ~ risk_score,
  data = tcga_final,
  reverse = TRUE
)$concordance

roc_tcga <- timeROC::timeROC(
  T = tcga_final$OS_time,
  delta = tcga_final$OS_event,
  marker = tcga_final$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(
    365,
    1095,
    1825
  )
)

c_nested <- survival::concordance(
  Surv(
    OS_time,
    OS_event
  ) ~ risk_score,
  data = nested_oof,
  reverse = TRUE
)$concordance

roc_nested <- timeROC::timeROC(
  T = nested_oof$OS_time,
  delta = nested_oof$OS_event,
  marker = nested_oof$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(
    365,
    1095,
    1825
  )
)

c_geo <- survival::concordance(
  Surv(
    OS_time_days,
    OS_event
  ) ~ risk_score,
  data = geo_validation,
  reverse = TRUE
)$concordance

roc_geo <- timeROC::timeROC(
  T = geo_validation$OS_time_days,
  delta = geo_validation$OS_event,
  marker = geo_validation$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(
    365,
    1095,
    1825
  )
)

performance_summary <- data.frame(

  Cohort = c(
    "TCGA final fit",
    "TCGA nested-CV",
    "GSE41613 external"
  ),

  N = c(
    nrow(tcga_final),
    nrow(nested_oof),
    nrow(geo_validation)
  ),

  Events = c(
    sum(tcga_final$OS_event),
    sum(nested_oof$OS_event),
    sum(geo_validation$OS_event)
  ),

  C_index = c(
    c_tcga,
    c_nested,
    c_geo
  ),

  AUC_1y = c(
    roc_tcga$AUC[1],
    roc_nested$AUC[1],
    roc_geo$AUC[1]
  ),

  AUC_3y = c(
    roc_tcga$AUC[2],
    roc_nested$AUC[2],
    roc_geo$AUC[2]
  ),

  AUC_5y = c(
    roc_tcga$AUC[3],
    roc_nested$AUC[3],
    roc_geo$AUC[3]
  )
)

write.table(
  performance_summary,
  "Final8Gene_Prognostic_Performance_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

auc_long <- bind_rows(

  data.frame(
    Cohort = performance_summary$Cohort,
    Time = "1-year",
    AUC = performance_summary$AUC_1y
  ),

  data.frame(
    Cohort = performance_summary$Cohort,
    Time = "3-year",
    AUC = performance_summary$AUC_3y
  ),

  data.frame(
    Cohort = performance_summary$Cohort,
    Time = "5-year",
    AUC = performance_summary$AUC_5y
  )
)

auc_long$Time <- factor(
  auc_long$Time,
  levels = c(
    "1-year",
    "3-year",
    "5-year"
  )
)

p_auc <- ggplot(
  auc_long,
  aes(
    x = Time,
    y = AUC,
    group = Cohort,
    shape = Cohort,
    linetype = Cohort
  )
) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted"
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        AUC
      )
    ),
    vjust = -0.8,
    size = 3.5,
    show.legend = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      0.45,
      0.80
    )
  ) +
  labs(
    title = "Time-dependent AUC across validation settings",
    x = NULL,
    y = "AUC",
    shape = NULL,
    linetype = NULL
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "bottom"
  )

ggsave(
  "Figure_Final8Gene_AUC_Comparison_PUBLICATION.pdf",
  p_auc,
  width = 9,
  height = 6
)

ggsave(
  "Figure_Final8Gene_AUC_Comparison_PUBLICATION.png",
  p_auc,
  width = 9,
  height = 6,
  dpi = 600
)

fit_risk <- survival::coxph(
  Surv(
    OS_time,
    OS_event
  ) ~ risk_z,
  data = clinical_dat
)

fit_clinical <- survival::coxph(
  Surv(
    OS_time,
    OS_event
  ) ~
    age +
    sex +
    stage4,
  data = clinical_dat
)

get_cindex <- function(fit, dat) {

  lp <- predict(
    fit,
    newdata = dat,
    type = "lp"
  )

  survival::concordance(
    Surv(
      dat$OS_time,
      dat$OS_event
    ) ~ lp,
    reverse = TRUE
  )$concordance
}

clinical_compare <- data.frame(

  Model = c(
    "8-gene risk",
    "Clinical",
    "Clinical + 8-gene risk"
  ),

  C_index = c(
    get_cindex(
      fit_risk,
      clinical_dat
    ),

    get_cindex(
      fit_clinical,
      clinical_dat
    ),

    get_cindex(
      combined_fit,
      clinical_dat
    )
  )
)

clinical_compare$Model <- factor(
  clinical_compare$Model,
  levels = c(
    "Clinical",
    "8-gene risk",
    "Clinical + 8-gene risk"
  )
)

p_cindex <- ggplot(
  clinical_compare,
  aes(
    x = Model,
    y = C_index
  )
) +
  geom_col(
    width = 0.62
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        C_index
      )
    ),
    vjust = -0.5,
    size = 4
  ) +
  coord_cartesian(
    ylim = c(
      0,
      0.78
    )
  ) +
  labs(
    title = "Incremental prognostic value of the 8-gene signature",
    x = NULL,
    y = "C-index"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 11
    )
  )

ggsave(
  "Figure_Final8Gene_Clinical_Cindex_PUBLICATION.pdf",
  p_cindex,
  width = 8,
  height = 6
)

ggsave(
  "Figure_Final8Gene_Clinical_Cindex_PUBLICATION.png",
  p_cindex,
  width = 8,
  height = 6,
  dpi = 600
)

cat("\n========================================\n")
cat("FINAL MODEL PERFORMANCE\n")
cat("========================================\n")

print(
  performance_summary
)

cat("\nClinical comparison:\n")

print(
  clinical_compare
)

cat("\nFrozen standardized formula:\n")

cat(
  formula_text,
  "\n"
)

cat("========================================\n")
