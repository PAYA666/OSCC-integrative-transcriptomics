# Analysis map

| Manuscript component | Main source location | Selected results |
|---|---|---|
| Original three-part analyses | historical/Part1, Part2, Part3 | Relevant frozen results below; original large project objects not bundled |
| Strict feature selection and internal validation | revision/02_strict_cv | results/02_strict_cv and results/selection_null |
| Smoking/HPV/alcohol/immune-marker sensitivities | revision/03_smoking | results/HPV_gene_level_sensitivity, HPV_score_sensitivity, Smoking_confounder_sensitivity, alcohol_matched, clinical_plus_MCPcounter |
| Resistance robustness | revision/04_resistance_loo | results/resistance_loo |
| External survival and treatment checks | revision/05_GSE41613_gene_level, 05_GSE42743, 09_GSE85446 | results/05_GSE41613_gene_level, 05_GSE42743 and Final8_* summary tables (GSE85446) |
| Co-expression conditioning | revision/11_network | results/all_score_pairs_by_module.tsv; all_module_hallmarks.tsv |
| Anatomy-restricted analyses | revision/13_site_de, 15_site_retrain | results/site_sensitivity; verified_sensitivities.tsv |
| Published 13-gene comparator | revision/14_benchmark | parameters/published_table2_parameters.tsv; results/GSE42743_paired_signature_comparison.tsv |
| SEQC technical-reference comparison | revision/21_R2_Major2_SEQC_bridge | results/R2_Major2_SEQC_bridge; results/S29_S30_machine_readable |
| GSE84846 deposited binary fields | revision/24_GSE84846 | results/S29_S30_machine_readable/Supplementary_Table_S29_GSE84846_binary_status.tsv |
| E-MTAB-15148 Ultra7 trajectories | historical/Part3/supp_step7_emtab15148_corrected.R | results/single_cell/Supplementary_Table_S31_EMTAB15148_Ultra7_trajectory.tsv |
| Figure generation | revision/06_figures | See manuscript figures; fix_figure5_title.py preserves data regions |

Consult result_aliases.tsv where duplicate frozen tables were retained under one canonical path. Data-file checksums record content correspondence; they do not certify execution of all source scripts.
