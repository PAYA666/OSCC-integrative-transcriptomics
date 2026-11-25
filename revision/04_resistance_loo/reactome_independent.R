suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
  library(msigdbr)
})

root <- Sys.getenv("OSCC_DATA", Sys.getenv("OSCC_DATA"))
base <- file.path(root, "Part3_Cisplatin", "Step1_FINAL_V4_OFFLINE")
sets <- msigdbr(species = "Homo sapiens", collection = "C2",
                subcollection = "CP:REACTOME")
pathways <- split(sets$gene_symbol, sets$gs_name)
pathways <- lapply(pathways, unique)
setDTthreads(1)

key_pathways <- c(
  "REACTOME_TNF_SIGNALING",
  "REACTOME_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY",
  "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",
  "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES"
)
not_tested_pathways <- "REACTOME_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR"
stopifnot(all(key_pathways %in% names(pathways)))

fits <- list()
eligibility <- list()
for (id in c("GSE197561", "GSE111585")) {
  fn <- function(s) file.path(base, id, paste0(id, s))
  x <- readRDS(fn("_GeneLevel_Expression.rds"))
  meta <- fread(fn("_PCA_coordinates.tsv"))[match(colnames(x), Sample)]
  stopifnot(identical(meta$Sample, colnames(x)),
            all(table(meta$Group) == 3L), !anyDuplicated(rownames(x)))
  original <- fread(fn("_limma_ALL.tsv"))
  eligibility[[length(eligibility) + 1L]] <- data.table(
    cohort = id, pathway = not_tested_pathways,
    overlap_size = sum(pathways[[not_tested_pathways]] %in% rownames(x)),
    status = "not_tested_below_minSize_15")

  for (drop_i in 0:6) {
    use <- if (drop_i == 0L) seq_len(ncol(x)) else setdiff(seq_len(ncol(x)), drop_i)
    group <- factor(meta$Group[use], levels = c("Parental", "Resistant"))
    fit <- eBayes(lmFit(x[, use, drop = FALSE], model.matrix(~group)))
    tab <- as.data.table(topTable(fit, coef = 2L, number = Inf,
                                  sort.by = "none"), keep.rownames = "Gene")
    if (drop_i == 0L) {
      check <- merge(tab, original, by = "Gene", suffixes = c(".new", ".old"))
      stopifnot(nrow(check) == nrow(original),
                max(abs(check$logFC.new - check$logFC.old)) < 1e-7,
                max(abs(check$t.new - check$t.old)) < 1e-6)
    }
    ord <- order(-tab$t, tab$Gene)
    ranks <- setNames(tab$t[ord], tab$Gene[ord])
    set.seed(20261122L + drop_i + if (id == "GSE111585") 100L else 0L)
    fg <- fgseaMultilevel(pathways = pathways, stats = ranks,
                          minSize = 15L, maxSize = 500L, eps = 0,
                          nproc = 1L)
    fg <- as.data.table(fg)[, .(pathway, size, ES, NES, pval, padj)]
    fg[, `:=`(cohort = id,
              fit = if (drop_i == 0L) "full" else paste0("drop_", meta$Sample[drop_i]),
              dropped_sample = if (drop_i == 0L) NA_character_ else meta$Sample[drop_i])]
    fits[[length(fits) + 1L]] <- fg
    fwrite(rbindlist(fits, use.names = TRUE),
           "reactome_all_fits.partial.tsv", sep = "\t")
  }
}

all_fits <- rbindlist(fits, use.names = TRUE)
stopifnot(uniqueN(all_fits[, .(cohort, fit)]) == 14L,
          all(key_pathways %in% all_fits$pathway))
key_summary <- all_fits[pathway %in% key_pathways][
  , .(fits = .N,
      positive_NES = sum(NES > 0),
      FDR05 = sum(padj < 0.05),
      nominal_p05 = sum(pval < 0.05),
      full_NES = NES[fit == "full"],
      full_padj = padj[fit == "full"]),
  by = .(cohort, pathway)]

fwrite(all_fits, "reactome_all_fits.tsv", sep = "\t")
fwrite(key_summary, "reactome_keypathway_stability.tsv", sep = "\t")
fwrite(rbindlist(eligibility), "reactome_pathway_eligibility.tsv", sep = "\t")
capture.output(sessionInfo(), file = "reactome_session.txt")
