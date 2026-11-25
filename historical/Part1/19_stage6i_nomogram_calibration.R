rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
setwd(project_root)

library(survival)
library(rms)
library(dplyr)

model_dat <- readRDS(
  "TCGA_OSCC_Nature8_Final8Gene_ClinicalCompleteCases.rds"
)

model_dat$sex <- factor(
  model_dat$sex,
  levels = c(
    "female",
    "male"
  )
)

model_dat$stage4 <- factor(
  model_dat$stage4,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IVA",
    "Stage IVB"
  )
)

stopifnot(
  "risk_z" %in% colnames(model_dat)
)

dd <- datadist(
  model_dat
)

options(
  datadist = "dd"
)

fit_nomogram <- rms::cph(
  Surv(
    OS_time,
    OS_event
  ) ~
    age +
    sex +
    stage4 +
    risk_z,

  data = model_dat,

  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 1095
)

surv_fun <- Survival(
  fit_nomogram
)

nom <- nomogram(
  fit_nomogram,

  fun = list(

    function(x)
      surv_fun(
        365,
        x
      ),

    function(x)
      surv_fun(
        1095,
        x
      ),

    function(x)
      surv_fun(
        1825,
        x
      )
  ),

  funlabel = c(
    "1-year OS",
    "3-year OS",
    "5-year OS"
  ),

  lp = FALSE,

  fun.at = c(
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9
  )
)

pdf(
  "Figure_Nature8_OSCC_Final8Gene_Nomogram.pdf",
  width = 13,
  height = 8
)

plot(
  nom,
  xfrac = 0.35,
  cex.axis = 0.8,
  cex.var = 0.9
)

dev.off()

png(
  "Figure_Nature8_OSCC_Final8Gene_Nomogram.png",
  width = 3200,
  height = 1900,
  res = 250
)

plot(
  nom,
  xfrac = 0.35,
  cex.axis = 0.8,
  cex.var = 0.9
)

dev.off()

set.seed(20260813)

cal_1y <- calibrate(
  fit_nomogram,
  cmethod = "KM",
  method = "boot",
  u = 365,
  m = 60,
  B = 500
)

set.seed(20260814)

cal_3y <- calibrate(
  fit_nomogram,
  cmethod = "KM",
  method = "boot",
  u = 1095,
  m = 60,
  B = 500
)

set.seed(20260815)

cal_5y <- calibrate(
  fit_nomogram,
  cmethod = "KM",
  method = "boot",
  u = 1825,
  m = 60,
  B = 500
)

pdf(
  "Figure_Nature8_OSCC_Final8Gene_Calibration.pdf",
  width = 13,
  height = 4.5
)

par(
  mfrow = c(1, 3),
  mar = c(5, 5, 4, 2)
)

plot(
  cal_1y,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "Predicted 1-year OS",
  ylab = "Observed 1-year OS",
  subtitles = FALSE
)

abline(
  0,
  1,
  lty = 2
)

title(
  "1-year calibration"
)

plot(
  cal_3y,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "Predicted 3-year OS",
  ylab = "Observed 3-year OS",
  subtitles = FALSE
)

abline(
  0,
  1,
  lty = 2
)

title(
  "3-year calibration"
)

plot(
  cal_5y,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "Predicted 5-year OS",
  ylab = "Observed 5-year OS",
  subtitles = FALSE
)

abline(
  0,
  1,
  lty = 2
)

title(
  "5-year calibration"
)

dev.off()

saveRDS(
  fit_nomogram,
  "TCGA_OSCC_Nature8_Final8Gene_NomogramModel.rds"
)

saveRDS(
  cal_1y,
  "TCGA_OSCC_Nature8_Final8Gene_Calibration_1y.rds"
)

saveRDS(
  cal_3y,
  "TCGA_OSCC_Nature8_Final8Gene_Calibration_3y.rds"
)

saveRDS(
  cal_5y,
  "TCGA_OSCC_Nature8_Final8Gene_Calibration_5y.rds"
)

cat("\n========================================\n")
cat("NOMOGRAM MODEL\n")
cat("========================================\n")

print(
  fit_nomogram
)

cat(
  "\nPatients:",
  nrow(model_dat),
  "\n"
)

cat(
  "Events:",
  sum(model_dat$OS_event),
  "\n"
)

cat("========================================\n")
