suppressPackageStartupMessages({library(DESeq2);library(data.table);library(fgsea)})
root <- Sys.getenv("OSCC_DATA", unset = NA_character_)
if (is.na(root) || !dir.exists(root)) stop("Set OSCC_DATA to the project data directory")
revision_root <- Sys.getenv("OSCC_REVISION", unset = normalizePath("..", winslash = "/", mustWork = TRUE))
dds0 <- readRDS(file.path(root,"Part2_TCGA_OSCC_Smoking_DESeq2_INPUT.rds"))
m <- as.data.table(as.data.frame(colData(dds0)))
man <- fread(file.path(revision_root,"01_baseline","results","sample_manifest.tsv"))[sample_type=="Primary Tumor"]
stopifnot(ncol(dds0)==332L,nrow(dds0)==17389L,!anyDuplicated(m$patient_id),!anyDuplicated(man$patient_barcode))
idx <- match(m$patient_id,man$patient_barcode)
stopifnot(!anyNA(idx),identical(as.character(m$smoking_bin),as.character(man$smoking_bin[idx])))
m[,site:=factor(fcase(grepl("Base of tongue",man$primary_site[idx]),"base_tongue",grepl("pharynx",man$primary_site[idx]),"mixed",grepl("tongue",man$primary_site[idx]),"other_tongue",grepl("Floor of mouth",man$primary_site[idx]),"floor",default="other_oral"))]
m[,sex:=factor(Sex)];m[,stage:=factor(Stage)];m[,age10:=as.numeric(as.character(Age))/10];m[,smoke:=factor(smoking_bin,levels=c("Never","Ever"))]
mcp <- readRDS(file.path(root,"Part2_Smoking_MCPcounter_Scores_OFFLINE.rds"))
stopifnot(is.matrix(mcp),nrow(mcp)==10L,!anyDuplicated(colnames(mcp)))
mcpdt <- as.data.table(t(mcp),keep.rownames="sample_barcode")
mcpdt[,patient_id:=substr(sample_barcode,1,12)]
stopifnot(!anyDuplicated(mcpdt$patient_id))
mcpdt[,sample_barcode:=NULL]
cellvars <- paste0("MCP_",make.names(rownames(mcp)))
setnames(mcpdt,old=rownames(mcp),new=cellvars)
m <- merge(m,mcpdt,by="patient_id",all.x=TRUE,sort=FALSE)
cellvars <- paste0("MCP_",make.names(rownames(mcp)))
basevars <- c("age10","sex","stage","site","smoke")
keep <- complete.cases(m[,c(basevars,cellvars),with=FALSE])
meta <- as.data.frame(droplevels(m[keep,c(basevars,cellvars),with=FALSE]))
ct <- counts(dds0)[,keep]
rownames(meta)<-colnames(ct)
stopifnot(identical(substr(colnames(ct),1,12),m$patient_id[keep]),!anyNA(meta),length(unique(meta$smoke))==2L);cat("complete N=",nrow(meta),"; ever=",sum(meta$smoke=="Ever"),"; never=",sum(meta$smoke=="Never"),"\\n",sep="")
meta[,cellvars] <- lapply(meta[,cellvars,drop=FALSE],function(z) as.numeric(scale(z)))
f2 <- ~age10+sex+stage+site+smoke
f3 <- reformulate(c("age10","sex","stage","site",cellvars,"smoke"))
X <- model.matrix(f3,meta)
cat("N",nrow(meta),"design_rank",qr(X)$rank,"columns",ncol(X),"condition",kappa(X),"\n")
stopifnot(qr(X)$rank==ncol(X))
fitone <- function(form,label){
 d <- DESeqDataSetFromMatrix(ct,meta,form)
 d <- DESeq(d,quiet=TRUE,parallel=FALSE)
 rr <- as.data.table(as.data.frame(results(d,contrast=c("smoke","Ever","Never"),alpha=.05)),keep.rownames="gene_id")
 main <- fread(file.path(root,"Part2_Smoking_DEG_Annotated.tsv"))
 rr[,symbol:=main$gene_name[match(gene_id,main$gene_id)]]
 rr[,model:=label]
 rr
}
clin <- fitone(f2,"clinical")
cell <- fitone(f3,"clinical_plus_MCPcounter")
stopifnot(identical(clin$gene_id,cell$gene_id))
fold <- merge(clin,cell,by="gene_id",suffixes=c("_clinical","_cell"))
common <- fold[padj_clinical<.05 & padj_cell<.05]
summary <- data.table(n=nrow(meta),never=sum(meta$smoke=="Never"),ever=sum(meta$smoke=="Ever"),
 cells_adjusted=paste(rownames(mcp),collapse=";"),genes=nrow(ct),
 clinical_FDR=sum(clin$padj<.05,na.rm=TRUE),cell_adjusted_FDR=sum(cell$padj<.05,na.rm=TRUE),
 both_FDR=nrow(common),same_direction=sum(sign(common$log2FoldChange_clinical)==sign(common$log2FoldChange_cell)),
 rho_all=cor(fold$log2FoldChange_clinical,fold$log2FoldChange_cell,method="spearman",use="complete.obs"),
 rho_common=if(nrow(common)>1)cor(common$log2FoldChange_clinical,common$log2FoldChange_cell,method="spearman") else NA_real_,
 design_rank=qr(X)$rank,design_columns=ncol(X),condition_number=kappa(X))
fwrite(clin,"clinical_reference.tsv",sep="\t");fwrite(cell,"clinical_plus_MCPcounter_genes.tsv",sep="\t");fwrite(summary,"clinical_plus_MCPcounter_summary.tsv",sep="\t")
cat("SUMMARY\n");print(summary)
gmt <- strsplit(readLines(file.path(root,"h.all.v2026.1.Hs.symbols.gmt")),"\t",fixed=TRUE)
pathways <- setNames(lapply(gmt,function(p)p[-(1:2)]),vapply(gmt,`[`,"",1L))
ranked <- function(z){bad<-duplicated(z$symbol)|duplicated(z$symbol,fromLast=TRUE);z<-z[!bad&!is.na(symbol)&nzchar(symbol)&is.finite(stat)];setorder(z,-stat,gene_id);setNames(z$stat,z$symbol)}
set.seed(20261123)
g <- fgseaMultilevel(pathways=pathways,stats=ranked(cell),minSize=15,maxSize=500,nproc=1,eps=0)
fwrite(g[,.(pathway,size,NES,pval,padj)],"clinical_plus_MCPcounter_hallmarks.tsv",sep="\t")
fwrite(data.table(input=c("Part2_TCGA_OSCC_Smoking_DESeq2_INPUT.rds","Part2_Smoking_MCPcounter_Scores_OFFLINE.rds","Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv",file.path(revision_root,"01_baseline","results","sample_manifest.tsv")),md5=unname(tools::md5sum(c(file.path(root,"Part2_TCGA_OSCC_Smoking_DESeq2_INPUT.rds"),file.path(root,"Part2_Smoking_MCPcounter_Scores_OFFLINE.rds"),file.path(root,"Part2_TCGA_OSCC_Nature8_Metadata_ExposureClinical.tsv"),file.path(revision_root,"01_baseline","results","sample_manifest.tsv"))))),"clinical_plus_MCPcounter_input_hashes.tsv",sep="\t")
capture.output(sessionInfo(),file="clinical_plus_MCPcounter_session.txt")

