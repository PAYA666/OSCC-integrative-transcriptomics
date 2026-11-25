root <- Sys.getenv("OSCC_DATA")
fig <- file.path(root,"Integration_Manuscript/Figure5_SharedNetworkArchitecture")
d <- read.delim(file.path(fig,"Figure5_MasterData.tsv"))
scores <- c("Stromal5","Proliferation2","Ultra7","Priority52","Exposure2","Final8Risk")
mes <- grep("^ME",names(d),value=TRUE)
stopifnot(!anyDuplicated(d$PatientID),length(mes)==7,nrow(d)==341)
out <- function(x,f) write.table(x,f,sep="\t",row.names=FALSE,quote=FALSE)
res <- do.call(rbind,lapply(combn(scores,2,simplify=FALSE),function(pair) {
  v <- na.omit(d[c(pair,mes)])
  x <- rank(v[[1]]); y <- rank(v[[2]])
  do.call(rbind,lapply(mes,function(me){
    z <- rank(v[[me]])
    data.frame(score1=pair[1],score2=pair[2],module=me,N=nrow(v),
      raw_rho=cor(x,y),rho_1_module=cor(x,z),rho_2_module=cor(y,z),
      partial_rho=cor(resid(lm(x~z)),resid(lm(y~z))))
  }))
}))
stopifnot(nrow(res)==105,all(is.finite(res$partial_rho)))
out(res,"all_score_pairs_by_module.tsv")
highlight <- subset(res,(score1=="Stromal5" & score2 %in% c("Ultra7","Priority52") & module=="MEturquoise") |
                  (score1=="Proliferation2" & score2=="Final8Risk" & module=="MEblack"))
print(highlight,row.names=FALSE)
src <- read.delim(file.path(root,"TCGA_OSCC_Nature8_WGCNA_gene_module_table.tsv"),na.strings=c("","NA"))
stopifnot(nrow(src)==5000,!anyDuplicated(src$gene_id))

ok <- !is.na(src$gene_name) & !duplicated(src$gene_name) & !duplicated(src$gene_name,fromLast=TRUE)
g <- src[ok,]
sets <- strsplit(readLines(file.path(root,"h.all.v2026.1.Hs.symbols.gmt")),"\t",fixed=TRUE)
modules <- setdiff(unique(g$module_color),"grey")
stopifnot(length(sets)==50,setequal(paste0("ME",modules),mes))
ora <- do.call(rbind,lapply(modules,function(mod){
  ids <- g$gene_name[g$module_color==mod]
  do.call(rbind,lapply(sets,function(set){
    path <- intersect(set[-(1:2)],g$gene_name)
    k <- length(intersect(ids,path)); N <- nrow(g)
    data.frame(module=mod,hallmark=set[1],overlap=k,module_n=length(ids),
               pathway_in_background=length(path),background=N,
               P=phyper(k-1,length(path),N-length(path),length(ids),lower.tail=FALSE))
  }))
}))
ora$FDR <- p.adjust(ora$P,method="BH")
ora <- ora[order(ora$module,ora$P),]
out(ora,"all_module_hallmarks.tsv")
print(ora[!duplicated(ora$module),c("module","hallmark","overlap","background","FDR")],row.names=FALSE)
