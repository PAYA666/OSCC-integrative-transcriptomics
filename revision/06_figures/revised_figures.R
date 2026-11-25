suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- ".."
out <- "."
fig_dir <- out

blue <- "#3C78B5"
red <- "#E05A5A"
amber <- "#E79B00"
slate <- "#788DA2"
dark <- "#30343B"
base_theme <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 10, color = dark),
    axis.text = element_text(size = 8.5, color = dark),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.position = "bottom",
    plot.margin = margin(4, 8, 4, 4)
  )

frequency <- fread(file.path(root, "02_strict_cv/results/strict_feature_frequency.tsv"))
final8 <- c("CTSG", "P4HA1", "NKAIN1", "LRRC38", "SPSB4", "SPINK6", "TNFRSF4", "DEFB1")
frequency <- frequency[gene %in% final8]
frequency[, gene := factor(gene, levels = rev(final8))]
p_a <- ggplot(frequency, aes(selection_frequency, gene)) +
  geom_vline(xintercept = 0.80, linetype = 2, linewidth = 0.45, color = "grey55") +
  geom_segment(aes(x = 0.5, xend = selection_frequency, yend = gene),
               linewidth = 0.8, color = slate) +
  geom_point(size = 2.6, color = blue) +
  geom_text(aes(label = sprintf("%.2f", selection_frequency)),
            hjust = -0.45, size = 2.7, color = dark) +
  scale_x_continuous(limits = c(0.5, 1.04), breaks = seq(0.5, 1.0, 0.1),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Selection frequency across 100 outer models", y = NULL) +
  base_theme + theme(legend.position = "none")

gene_hr <- fread(file.path(root, "05_GSE41613_gene_level/GSE41613_individual_gene_HR.tsv"))
gene_hr[, gene := factor(gene, levels = rev(final8))]
gene_hr[, direction := ifelse(HR_per_external_SD > 1, "Higher hazard", "Lower hazard")]
gene_hr[, label := sprintf("%.2f (%.2f–%.2f)", HR_per_external_SD, lower95, upper95)]
p_b <- ggplot(gene_hr, aes(HR_per_external_SD, gene, color = direction)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, color = "grey55") +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), width = 0.18,
                orientation = "y", linewidth = 0.65) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("Higher hazard" = red, "Lower hazard" = blue)) +
  scale_x_log10(limits = c(0.45, 3.8), breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  labs(x = "Hazard ratio per external-cohort SD", y = NULL) +
  base_theme

repeat_c <- fread(file.path(root, "02_strict_cv/results/repeat_cindex.tsv"))
repeat_c[, repeat_id := factor(repeat_id)]
p_c <- ggplot(repeat_c, aes(x = 1, y = C)) +
  geom_boxplot(width = 0.22, outlier.shape = NA, fill = "grey92",
               color = slate, linewidth = 0.55) +
  geom_jitter(width = 0.07, height = 0, size = 2.2, color = blue,
              shape = 21, fill = "white", stroke = 0.7) +
  geom_hline(yintercept = mean(repeat_c$C), color = amber, linewidth = 0.8) +
  annotate("text", x = 1.08, y = 0.612,
           label = sprintf("Mean %.3f (SD %.3f)", mean(repeat_c$C), sd(repeat_c$C)),
           hjust = 0, vjust = 1, size = 2.7, color = dark) +
  scale_x_continuous(breaks = 1, labels = "10 repeats") +
  coord_cartesian(ylim = c(0.54, 0.615), xlim = c(0.7, 1.55)) +
  labs(x = NULL, y = "Harrell's C-index") +
  base_theme + theme(legend.position = "none")

strict_c <- fread(file.path(root, "02_strict_cv/results/strict_pipeline_cindex.tsv"))
clinical_c <- fread(file.path(root, "02_strict_cv/results/clinical_vs_Final8_nested_cv.tsv"))
gse41613_c <- fread(file.path(root, "01_baseline/results/GSE41613_cindex_bootstrap.tsv"))
gse42743_c <- fread(file.path(root, "05_GSE42743/results/GSE42743_external_validation.tsv"))
c_data <- data.table(
  cohort = c("Strict genomic pipeline", "Clinical-only (same cases)",
             "Genomic pipeline (same cases)", "GSE41613, frozen Final8", "GSE42743, frozen Final8"),
  C = c(strict_c$C, clinical_c$clinical_C, clinical_c$Final8_C,
        gse41613_c$C_index, gse42743_c$C_index[1]),
  lower = c(strict_c$lower95_patient_bootstrap,
            clinical_c$clinical_C_lower95_patient_bootstrap,
            clinical_c$Final8_C_lower95_patient_bootstrap,
            gse41613_c$lower95_patient_bootstrap,
            gse42743_c$C_lower95_bootstrap[1]),
  upper = c(strict_c$upper95_patient_bootstrap,
            clinical_c$clinical_C_upper95_patient_bootstrap,
            clinical_c$Final8_C_upper95_patient_bootstrap,
            gse41613_c$upper95_patient_bootstrap,
            gse42743_c$C_upper95_bootstrap[1])
)
c_data[, cohort := factor(cohort, levels = rev(cohort))]
c_data[, label := sprintf("%.3f (%.3f–%.3f)", C, lower, upper)]
c_data[, group := ifelse(grepl("GSE", cohort), "External Final8", "Internal evaluation")]
p_d <- ggplot(c_data, aes(C, cohort, color = group)) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.45, color = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2,
                orientation = "y", linewidth = 0.65) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("Internal evaluation" = amber, "External Final8" = red)) +
  scale_x_continuous(limits = c(0.48, 0.86), breaks = seq(0.5, 0.8, 0.1)) +
  labs(x = "Harrell's C-index (95% patient-bootstrap CI)", y = NULL) +
  base_theme

strict_auc <- fread(file.path(root, "02_strict_cv/results/strict_pipeline_timeAUC.tsv"))
clinical_auc <- fread(file.path(root, "02_strict_cv/results/clinical_vs_Final8_timeAUC.tsv"))
gse42743_auc <- fread(file.path(root, "05_GSE42743/results/GSE42743_timeAUC.tsv"))
auc_data <- rbindlist(list(
  data.table(cohort = "Strict genomic pipeline", years = strict_auc$years,
             AUC = strict_auc$AUC, lower = strict_auc$lower95, upper = strict_auc$upper95),
  data.table(cohort = "GSE41613, frozen Final8", years = c(1, 3, 5),
             AUC = c(0.740, 0.707, 0.711), lower = c(0.631, 0.603, 0.602),
             upper = c(0.848, 0.811, 0.821)),
  data.table(cohort = "GSE42743, frozen Final8", years = gse42743_auc$years,
             AUC = gse42743_auc$AUC, lower = gse42743_auc$lower95,
             upper = gse42743_auc$upper95)
))
auc_data[, cohort := factor(cohort, levels = c("Strict genomic pipeline",
  "GSE41613, frozen Final8", "GSE42743, frozen Final8"))]
auc_data[, years := factor(years, levels = c(1, 3, 5), labels = c("1 year", "3 years", "5 years"))]
p_e <- ggplot(auc_data, aes(years, AUC, color = cohort, group = cohort, shape = cohort)) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.45, color = "grey65") +
  geom_line(linewidth = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Strict genomic pipeline" = amber,
                                "GSE41613, frozen Final8" = red,
                                "GSE42743, frozen Final8" = blue)) +
  scale_y_continuous(limits = c(0.35, 1), breaks = seq(0.4, 1, 0.2)) +
  labs(x = NULL, y = "Time-dependent AUC (95% CI)") + base_theme

gse41613_adj <- fread(file.path(root, "01_baseline/results/GSE41613_continuous_score.tsv"))
gse42743_adj <- fread(file.path(root, "05_GSE42743/results/GSE42743_adjusted_validation.tsv"))
adj_data <- data.table(
  cohort = c("GSE41613, adjusted", "GSE42743, adjusted"),
  HR = c(gse41613_adj$HR_per_external_SD[2], gse42743_adj$HR_per_SD),
  lower = c(gse41613_adj$lower95[2], gse42743_adj$lower95_HR),
  upper = c(gse41613_adj$upper95[2], gse42743_adj$upper95_HR),
  p = c(gse41613_adj$p_LR[2], gse42743_adj$p_LR_added_score)
)
adj_data[, cohort := factor(cohort, levels = rev(cohort))]
adj_data[, label := sprintf("HR %.2f (%.2f–%.2f), P=%.3f", HR, lower, upper, p)]
p_f <- ggplot(adj_data, aes(HR, cohort, color = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, color = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.18,
                orientation = "y", linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("GSE41613, adjusted" = red, "GSE42743, adjusted" = blue)) +
  scale_x_log10(limits = c(0.65, 3.9), breaks = c(0.75, 1, 1.5, 2, 3)) +
  labs(x = "Adjusted hazard ratio per cohort SD", y = NULL) +
  base_theme + theme(legend.position = "none")

p_b <- p_b + theme(legend.position = "none")
p_d <- p_d + theme(legend.position = "none")
p_e <- p_e + theme(legend.position = "bottom")

fig6 <- (p_a | p_b) / (p_c | p_d) / (p_e | p_f) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14, color = "black"))

comparison_c <- data.table(
  model = factor(c("Clinical-only", "Genomic pipeline"),
                 levels = c("Clinical-only", "Genomic pipeline")),
  C = c(clinical_c$clinical_C, clinical_c$Final8_C),
  lower = c(clinical_c$clinical_C_lower95_patient_bootstrap,
            clinical_c$Final8_C_lower95_patient_bootstrap),
  upper = c(clinical_c$clinical_C_upper95_patient_bootstrap,
            clinical_c$Final8_C_upper95_patient_bootstrap)
)
p7a <- ggplot(comparison_c, aes(C, model, color = model, shape = model)) +
  geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.45, color = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2,
                orientation = "y", linewidth = 0.75) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Clinical-only" = slate, "Genomic pipeline" = amber)) +
  scale_x_continuous(limits = c(0.48, 0.74), breaks = seq(0.5, 0.7, 0.05)) +
  labs(x = "Harrell's C-index (95% patient-bootstrap CI)", y = NULL) +
  base_theme + theme(legend.position = "none")

clinical_auc[, years := factor(years, levels = c(1, 3, 5), labels = c("1 year", "3 years", "5 years"))]
pipeline_auc <- clinical_auc[model == "score"]
clinical_auc_plot <- clinical_auc[model == "clinical_lp"]
p7b_data <- rbindlist(list(
  data.table(model = "Clinical-only", years = clinical_auc_plot$years,
             AUC = clinical_auc_plot$AUC, lower = clinical_auc_plot$lower95,
             upper = clinical_auc_plot$upper95),
  data.table(model = "Genomic pipeline", years = pipeline_auc$years,
             AUC = pipeline_auc$AUC, lower = pipeline_auc$lower95,
             upper = pipeline_auc$upper95)
))
p7b <- ggplot(p7b_data, aes(years, AUC, color = model, group = model, shape = model)) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.45, color = "grey65") +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08, linewidth = 0.6) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("Clinical-only" = slate, "Genomic pipeline" = amber)) +
  scale_y_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, 0.1)) +
  labs(x = NULL, y = NULL) +
  base_theme + theme(legend.position = "bottom")

fig7 <- (p7a | p7b) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14, color = "black"))

ggsave(file.path(fig_dir, "Fig6_revised.png"), fig6, width = 12, height = 9,
       units = "in", dpi = 300, bg = "white", device = grDevices::png,
       type = "cairo")
ggsave(file.path(fig_dir, "Fig6_revised.pdf"), fig6, width = 12, height = 9,
       units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(fig_dir, "Fig7_revised.png"), fig7, width = 11.6, height = 5.4,
       units = "in", dpi = 300, bg = "white", device = grDevices::png,
       type = "cairo")
ggsave(file.path(fig_dir, "Fig7_revised.pdf"), fig7, width = 11.6, height = 5.4,
       units = "in", device = cairo_pdf, bg = "white")
capture.output(sessionInfo(), file = file.path(out, "session.txt"))
