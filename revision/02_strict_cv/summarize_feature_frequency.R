suppressPackageStartupMessages(library(data.table))

paths <- list.files("checkpoints", pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(paths) == 100L)
fits <- lapply(paths, readRDS)
stopifnot(length(unique(vapply(fits, `[[`, integer(1), "repeat_id"))) == 10L)

selected <- rbindlist(lapply(fits, function(fit) {
  genes <- unname(fit$beta_symbols)
  data.table(repeat_id = fit$repeat_id, fold_id = fit$fold_id,
             gene = genes)
}))
frequency <- selected[, .(selected_outer_models = uniqueN(paste(repeat_id, fold_id)),
                           selection_frequency = uniqueN(paste(repeat_id, fold_id)) / 100),
                       by = "gene"][order(-selection_frequency, gene)]
fwrite(frequency, "results/strict_feature_frequency.tsv", sep = "\t")
print(frequency)
