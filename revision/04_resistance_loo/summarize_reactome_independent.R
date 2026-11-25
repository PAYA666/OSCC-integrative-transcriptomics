suppressPackageStartupMessages(library(data.table))

x <- fread("reactome_all_fits.tsv")
full <- x[fit == "full"]
a <- full[cohort == "GSE197561",
          .(pathway, NES_197561 = NES, pval_197561 = pval,
            padj_197561 = padj)]
b <- full[cohort == "GSE111585",
          .(pathway, NES_111585 = NES, pval_111585 = pval,
            padj_111585 = padj)]
cross <- merge(a, b, by = "pathway", all = FALSE)
cross[, same_direction := sign(NES_197561) == sign(NES_111585)]
cross[, both_FDR05 := padj_197561 < 0.05 & padj_111585 < 0.05]
cross[, worst_q := pmax(padj_197561, padj_111585)]
setorder(cross, -both_FDR05, -same_direction, worst_q)

fwrite(cross, "reactome_crosscohort_full.tsv", sep = "\t")
fwrite(cross[both_FDR05 & same_direction],
       "reactome_shared_FDR05_same_direction.tsv", sep = "\t")

stopifnot(uniqueN(x[, .(cohort, fit)]) == 14L,
          all(cross$pathway %in% full$pathway))
