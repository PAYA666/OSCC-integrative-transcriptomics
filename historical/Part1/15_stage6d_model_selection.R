rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(dplyr)

set.seed(20260813)

stability_table <- readRDS(
  "TCGA_OSCC_Nature8_RepeatedLASSO_Stability.rds"
)

datExpr <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_datExpr_top5000MAD.rds"
)

clinical_wgcna <- readRDS(
  "TCGA_OSCC_Nature8_WGCNA_clinical_aligned.rds"
)

keep <- (
  !is.na(clinical_wgcna$OS_time) &
  clinical_wgcna$OS_time > 0 &
  !is.na(clinical_wgcna$OS_event)
)

clinical_surv <- clinical_wgcna[
  keep, , drop = FALSE
]

expr_surv <- datExpr[
  keep, , drop = FALSE
]

genes4 <- stability_table$gene_id[
  stability_table$selection_frequency >= 0.80
]

genes5 <- stability_table$gene_id[
  stability_table$selection_frequency >= 0.60
]

genes14 <- stability_table$gene_id[
  stability_table$selection_frequency >= 0.50
]

cat("\nModel gene numbers:\n")

print(
  c(
    Model4 = length(genes4),
    Model5 = length(genes5),
    Model14 = length(genes14)
  )
)

gene_anno <- stability_table %>%
  dplyr::select(
    gene_id,
    gene_name,
    selection_frequency,
    module_color,
    MM,
    log2FoldChange,
    DEG_status
  )

uni_results <- list()

for (gene in genes14) {

  tmp <- data.frame(
    OS_time = clinical_surv$OS_time,
    OS_event = clinical_surv$OS_event,
    expr = expr_surv[, gene]
  )

  fit <- survival::coxph(
    survival::Surv(
      OS_time,
      OS_event
    ) ~ expr,
    data = tmp
  )

  s <- summary(fit)

  uni_results[[gene]] <- data.frame(
    gene_id = gene,
    beta = s$coefficients[1, "coef"],
    HR = s$coefficients[1, "exp(coef)"],
    lower95 = s$conf.int[1, "lower .95"],
    upper95 = s$conf.int[1, "upper .95"],
    P = s$coefficients[1, "Pr(>|z|)"]
  )
}

uni_cox <- dplyr::bind_rows(
  uni_results
)

uni_cox$FDR <- p.adjust(
  uni_cox$P,
  method = "BH"
)

uni_cox <- dplyr::left_join(
  uni_cox,
  gene_anno,
  by = "gene_id"
)

uni_cox <- uni_cox %>%
  dplyr::arrange(P)

make_model_data <- function(genes) {

  x <- as.data.frame(
    expr_surv[, genes, drop = FALSE]
  )

  colnames(x) <- genes

  data.frame(
    OS_time = clinical_surv$OS_time,
    OS_event = clinical_surv$OS_event,
    x,
    check.names = FALSE
  )
}

fit_cox_model <- function(genes) {

  dat <- make_model_data(genes)

  rhs <- paste(
    paste0("`", genes, "`"),
    collapse = " + "
  )

  f <- as.formula(
    paste(
      "Surv(OS_time, OS_event) ~",
      rhs
    )
  )

  survival::coxph(
    f,
    data = dat,
    x = TRUE,
    y = TRUE
  )
}

fit4 <- fit_cox_model(genes4)
fit5 <- fit_cox_model(genes5)
fit14 <- fit_cox_model(genes14)

extract_cox <- function(fit, model_name) {

  s <- summary(fit)

  out <- data.frame(
    gene_id = rownames(s$coefficients),
    beta = s$coefficients[, "coef"],
    HR = s$coefficients[, "exp(coef)"],
    lower95 = s$conf.int[, "lower .95"],
    upper95 = s$conf.int[, "upper .95"],
    P = s$coefficients[, "Pr(>|z|)"],
    Model = model_name,
    row.names = NULL
  )

  out
}

multi4 <- extract_cox(
  fit4,
  "4-gene"
)

multi5 <- extract_cox(
  fit5,
  "5-gene"
)

multi14 <- extract_cox(
  fit14,
  "14-gene"
)

multi4 <- dplyr::left_join(
  multi4,
  gene_anno[, c("gene_id","gene_name")],
  by = "gene_id"
)

multi5 <- dplyr::left_join(
  multi5,
  gene_anno[, c("gene_id","gene_name")],
  by = "gene_id"
)

multi14 <- dplyr::left_join(
  multi14,
  gene_anno[, c("gene_id","gene_name")],
  by = "gene_id"
)

get_cindex <- function(fit) {

  lp <- predict(
    fit,
    type = "lp"
  )

  survival::concordance(
    survival::Surv(
      fit$y[, 1],
      fit$y[, 2]
    ) ~ lp,
    reverse = TRUE
  )$concordance
}

c4 <- get_cindex(fit4)
c5 <- get_cindex(fit5)
c14 <- get_cindex(fit14)

bootstrap_corrected_cindex <- function(
    genes,
    B = 500,
    seed = 20260813) {

  set.seed(seed)

  dat <- make_model_data(genes)

  rhs <- paste(
    paste0("`", genes, "`"),
    collapse = " + "
  )

  form <- as.formula(
    paste(
      "Surv(OS_time, OS_event) ~",
      rhs
    )
  )

  fit_original <- survival::coxph(
    form,
    data = dat,
    x = TRUE,
    y = TRUE
  )

  lp_original <- predict(
    fit_original,
    newdata = dat,
    type = "lp"
  )

  C_original <- survival::concordance(
    survival::Surv(
      dat$OS_time,
      dat$OS_event
    ) ~ lp_original,
    reverse = TRUE
  )$concordance

  optimism <- rep(
    NA_real_,
    B
  )

  for (b in seq_len(B)) {

    idx <- sample(
      seq_len(nrow(dat)),
      replace = TRUE
    )

    boot_dat <- dat[
      idx,
      ,
      drop = FALSE
    ]

    fit_boot <- try(
      survival::coxph(
        form,
        data = boot_dat
      ),
      silent = TRUE
    )

    if (inherits(
      fit_boot,
      "try-error"
    )) {
      next
    }

    lp_boot <- predict(
      fit_boot,
      newdata = boot_dat,
      type = "lp"
    )

    C_boot <- survival::concordance(
      survival::Surv(
        boot_dat$OS_time,
        boot_dat$OS_event
      ) ~ lp_boot,
      reverse = TRUE
    )$concordance

    lp_test <- predict(
      fit_boot,
      newdata = dat,
      type = "lp"
    )

    C_test <- survival::concordance(
      survival::Surv(
        dat$OS_time,
        dat$OS_event
      ) ~ lp_test,
      reverse = TRUE
    )$concordance

    optimism[b] <- C_boot - C_test

    if (b %% 100 == 0) {

      cat(
        "Bootstrap:",
        b,
        "/",
        B,
        "\n"
      )
    }
  }

  mean_optimism <- mean(
    optimism,
    na.rm = TRUE
  )

  corrected_C <- C_original -
    mean_optimism

  data.frame(
    Apparent_Cindex = C_original,
    Mean_optimism = mean_optimism,
    Corrected_Cindex = corrected_C,
    Successful_bootstraps =
      sum(!is.na(optimism))
  )
}

cat("\nBootstrap Model 4...\n")

boot4 <- bootstrap_corrected_cindex(
  genes4,
  B = 500,
  seed = 20260814
)

cat("\nBootstrap Model 5...\n")

boot5 <- bootstrap_corrected_cindex(
  genes5,
  B = 500,
  seed = 20260815
)

cat("\nBootstrap Model 14...\n")

boot14 <- bootstrap_corrected_cindex(
  genes14,
  B = 500,
  seed = 20260816
)

model_compare <- data.frame(

  Model = c(
    "4-gene",
    "5-gene",
    "14-gene"
  ),

  Genes = c(
    length(genes4),
    length(genes5),
    length(genes14)
  ),

  AIC = c(
    AIC(fit4),
    AIC(fit5),
    AIC(fit14)
  ),

  Apparent_Cindex = c(
    c4,
    c5,
    c14
  ),

  Bootstrap_Cindex = c(
    boot4$Corrected_Cindex,
    boot5$Corrected_Cindex,
    boot14$Corrected_Cindex
  ),

  Optimism = c(
    boot4$Mean_optimism,
    boot5$Mean_optimism,
    boot14$Mean_optimism
  )
)

write.table(
  uni_cox,
  "TCGA_OSCC_Nature8_StableGenes_UnivariateCox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  multi4,
  "TCGA_OSCC_Nature8_Model4_MultivariableCox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  multi5,
  "TCGA_OSCC_Nature8_Model5_MultivariableCox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  multi14,
  "TCGA_OSCC_Nature8_Model14_MultivariableCox.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  model_compare,
  "TCGA_OSCC_Nature8_ModelComparison.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  fit4,
  "TCGA_OSCC_Nature8_Model4_Cox.rds"
)

saveRDS(
  fit5,
  "TCGA_OSCC_Nature8_Model5_Cox.rds"
)

saveRDS(
  fit14,
  "TCGA_OSCC_Nature8_Model14_Cox.rds"
)

saveRDS(
  model_compare,
  "TCGA_OSCC_Nature8_ModelComparison.rds"
)

cat("\n========================================\n")
cat("MODEL COMPARISON\n")
cat("========================================\n")

print(model_compare)
