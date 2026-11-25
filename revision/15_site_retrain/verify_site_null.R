suppressPackageStartupMessages(library(data.table))
out <- "results/site_selection_null_permutations"
f <- sort(list.files(out, pattern="^perm_[0-9]{4}\\.tsv$", full.names=TRUE))
stopifnot(length(f)==100L, identical(basename(f),sprintf("perm_%04d.tsv",1:100)))
genes <- c("CTSG","P4HA1","NKAIN1","LRRC38","SPSB4","SPINK6","TNFRSF4","DEFB1")
a <- rbindlist(lapply(f,fread))
stopifnot(nrow(a)==800L, setequal(a$gene,genes))
stopifnot(all(a$selected_outer_models %in% 0:100),
          all(abs(a$selection_frequency-a$selected_outer_models/100)<1e-12))
id <- as.integer(sub(".*perm_([0-9]{4})\\.tsv$","\\1",f))
for(i in seq_along(f)) stopifnot(all(a[permutation==id[i],.N]==8L), all(a[permutation==id[i],gene] %chin% genes))
s <- fread("results/site_selection_null_summary.tsv")
o <- fread("results/site_selection_null_overall.tsv")
screen <- fread("results/selection_screen_comparison.tsv")
stopifnot(nrow(s)==8L,nrow(o)==1L,o$permutations==100L,nrow(screen)==100L)
any <- a[,.(hit=any(selection_frequency>=.8)),by=permutation]
stopifnot(o$permutations_with_any_final8_gene_ge_0_80==sum(any$hit),
          abs(o$null_probability_any_final8_gene_ge_0_80-mean(any$hit))<1e-12)
for(g in genes){
  x <- a[gene==g][order(selection_frequency)]
  y <- s[gene==g]
  stopifnot(nrow(x)==100L,nrow(y)==1L,
    abs(y$median_selection_frequency-median(x$selection_frequency))<1e-12,
    abs(y$q95_selection_frequency-unname(quantile(x$selection_frequency,.95)))<1e-12,
    y$permutations_ge_0_80==sum(x$selection_frequency>=.8))
}
ck <- list.files("checkpoints",pattern="^r[0-9]{2}f[0-9]{2}\\.rds$",full.names=TRUE)
stopifnot(length(ck)==100L)
observed <- rbindlist(lapply(ck,function(p){
  z<-readRDS(p); data.table(repeat_id=z$repeat_id,fold_id=z$fold_id,gene=unname(z$beta_symbols))
}))[,.(selected_outer_models=.N,selection_frequency=.N/100),by=gene]
observed <- merge(data.table(gene=genes),observed,by="gene",all.x=TRUE)
observed[is.na(selected_outer_models),`:=`(selected_outer_models=0L,selection_frequency=0)]
fwrite(observed,"results/site_final8_observed_selection.tsv",sep="\t")
print(list(permutations=nrow(f),genes=length(genes),
           any_gene_ge_0_80=o$permutations_with_any_final8_gene_ge_0_80,
           null_probability=o$null_probability_any_final8_gene_ge_0_80,
           score_wald_agreement=mean(screen$threshold_agreement),observed=observed))
