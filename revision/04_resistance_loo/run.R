library(data.table)
library(limma)
library(fgsea)
root <- Sys.getenv("OSCC_DATA",Sys.getenv("OSCC_DATA"))
base <- file.path(root,"Part3_Cisplatin","Step1_FINAL_V4_OFFLINE")
gmt <- gmtPathways(file.path(root,"h.all.v2026.1.Hs.symbols.gmt"))
stopifnot(length(gmt)==50L)
ultra <- c("CCN1","CCL2","ZBTB10","DAB2","JUN","KYNU","CCL20")
setDTthreads(1)
genes <- paths <- checks <- list()
for(id in c("GSE197561","GSE111585")){
 fn <- function(s) file.path(base,id,paste0(id,s))
 x <- readRDS(fn("_GeneLevel_Expression.rds"))
 m <- fread(fn("_PCA_coordinates.tsv"))[match(colnames(x),Sample)]
 orig <- fread(fn("_limma_ALL.tsv"))
 stopifnot(identical(m$Sample,colnames(x)),all(table(m$Group)==3),
           !anyDuplicated(rownames(x)),all(ultra %in% rownames(x)))
 for(i in 0:6){
   use <- setdiff(1:6,i)
   g <- factor(m$Group[use],levels=c("Parental","Resistant"))
   fit <- eBayes(lmFit(x[,use],model.matrix(~g)))
   tab <- as.data.table(topTable(fit,coef=2,number=Inf,sort.by="none"),keep.rownames="Gene")
   if(!i){
     comp <- merge(tab,orig,by="Gene")
     stopifnot(nrow(comp)==nrow(orig),max(abs(comp$logFC.x-comp$logFC.y))<1e-7,
               max(abs(comp$t.x-comp$t.y))<1e-6)
     full <- tab
   }
   drop <- if(i) m$Sample[i] else "none"
   z <- merge(full[,.(Gene,full_lfc=logFC,full_q=adj.P.Val)],tab,by="Gene")
   z[,strict:=adj.P.Val<.05 & abs(logFC)>1]
   q <- z[full_q<.05 & abs(full_lfc)>1]
   checks[[length(checks)+1]] <- data.table(id,drop,full_strict=nrow(q),
     retained_direction=sum(sign(q$logFC)==sign(q$full_lfc)),retained_strict_same=sum(q$strict & sign(q$logFC)==sign(q$full_lfc)),
     genomewide_rho=cor(z$logFC,z$full_lfc,method="spearman"))
   genes[[length(genes)+1]] <- z[Gene %in% ultra,.(id,drop,Gene,logFC,adj.P.Val,strict)]
   ranks <- sort(setNames(tab$t,tab$Gene),decreasing=TRUE)
   set.seed(20261120+i)
   fg <- fgseaSimple(gmt,ranks,nperm=10000,minSize=15,maxSize=500,nproc=1)
   paths[[length(paths)+1]] <- fg[,.(id,drop,pathway,NES,pval,padj)]
 }
}
g <- rbindlist(genes);p <- rbindlist(paths);q <- rbindlist(checks)
stopifnot(nrow(g)==98, nrow(q)==14, nrow(p)==700)
fwrite(g,"ultra7_all_fits.tsv",sep="\t")
fwrite(p,"hallmark_all_fits.tsv",sep="\t")
fwrite(q,"original_DEG_retention.tsv",sep="\t")
fwrite(g[drop!="none",.(n=.N,min_logFC=min(logFC),max_logFC=max(logFC),positive=sum(logFC>0),strict=sum(strict)),by=.(id,Gene)],"ultra7_summary.tsv",sep="\t")
fwrite(p[drop!="none",.(n=.N,min_NES=min(NES),max_NES=max(NES),positive=sum(NES>0),FDR05=sum(padj<.05)),by=.(id,pathway)],"hallmark_summary.tsv",sep="\t")
capture.output(sessionInfo(),file="session.txt")
