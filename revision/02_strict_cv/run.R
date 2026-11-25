suppressPackageStartupMessages({
 library(DESeq2); library(WGCNA); library(glmnet); library(survival); library(data.table)
})
root <- Sys.getenv("OSCC_DATA",Sys.getenv("OSCC_DATA"))
b <- "../01_baseline/results"
meta <- fread(file.path(b,"sample_manifest.tsv"),na.strings="NA")
folds <- fread(file.path(b,"outer_folds.tsv"))
v <- meta[sample_type=="Primary Tumor" & !is.na(OS_event)]
stopifnot(nrow(v)==297L,!anyDuplicated(v$patient_barcode))
raw <- readRDS(file.path(root,"TCGA_OSCC_Nature8_raw_count_matrix.rds"))
raw <- raw[raw$gene_type=="protein_coding",]
cnt <- as.matrix(raw[,meta$sample_barcode])
rownames(cnt) <- raw$gene_id
sym <- setNames(raw$gene_name,raw$gene_id)
stopifnot(ncol(cnt)==373L,nrow(cnt)==19962L,!anyDuplicated(rownames(cnt)),all(cnt==round(cnt)))
storage.mode(cnt) <- "integer"
rm(raw)
args <- as.integer(commandArgs(TRUE))
stopifnot(length(args)==2,1<=args[1],args[1]<=args[2],args[2]<=100)
dir.create("checkpoints",showWarnings=FALSE)
setDTthreads(2)
logmsg <- function(...) {cat(format(Sys.time()),..., "\n");flush.console()}
stamp <- unname(tools::md5sum(c("run.R",file.path(b,c("outer_folds.tsv","sample_manifest.tsv")),
                                 file.path(root,"TCGA_OSCC_Nature8_raw_count_matrix.rds"))))
fitfold <- function(r,k){
 seed <- 20261120L+100L*r+k
 set.seed(seed)
 hold <- folds[repeat_id==r & fold_id==k,patient_barcode]
 tr <- meta[!patient_barcode %in% hold]
 te <- v[match(hold,patient_barcode)]
 va <- v[!patient_barcode %in% hold]
 stopifnot(!any(hold %in% tr$patient_barcode),!anyNA(te$patient_barcode),all(va$sample_barcode %in% tr$sample_barcode))
 c0 <- cnt[,tr$sample_barcode,drop=FALSE]
 c0 <- c0[rowSums(c0>=10L)>=10L,,drop=FALSE]
 typ <- factor(ifelse(tr$sample_type=="Primary Tumor","Tumor","Normal"),levels=c("Normal","Tumor"))
 dd <- DESeqDataSetFromMatrix(c0,data.frame(type=typ,row.names=colnames(c0)),~type)
 logmsg("DESeq2",r,k)
 dd <- DESeq(dd,quiet=TRUE,parallel=FALSE)
 rr <- results(dd,contrast=c("type","Tumor","Normal"),alpha=0.05)
 deg <- rownames(rr)[which(rr$padj<0.05 & abs(rr$log2FoldChange)>1)]
 fun <- dispersionFunction(dd)
 stopifnot(attr(fun,"fitType")=="parametric")
 ab <- attr(fun,"coefficients")
 trans <- function(q) log2((1+ab["extraPois"]+2*ab["asymptDisp"]*q+
                             2*sqrt(ab["asymptDisp"]*q*(1+ab["extraPois"]+ab["asymptDisp"]*q)))/(4*ab["asymptDisp"]))
 lg <- rowMeans(log(c0))
 sf <- function(z) apply(z,2,function(col) {
   keep <- is.finite(lg) & col>0
   exp(median(log(col[keep])-lg[keep]))
 })
 stopifnot(max(abs(sf(c0)-sizeFactors(dd))) < 1e-6)
 normtr <- trans(counts(dd,normalized=TRUE))
 chk <- assay(varianceStabilizingTransformation(dd,blind=FALSE))
 stopifnot(max(abs(normtr-chk))<1e-7)
 rm(chk,dd,rr)
 xt <- t(normtr[,typ=="Tumor",drop=FALSE])
 mads <- matrixStats::colMads(xt)
 top <- colnames(xt)[order(-mads,colnames(xt))[seq_len(5000)]]
 logmsg("WGCNA",r,k)
 net <- blockwiseModules(xt[,top],power=8,networkType="signed",TOMType="signed",corType="bicor",
                         minModuleSize=30,deepSplit=2,mergeCutHeight=.25,pamRespectsDendro=FALSE,
                         numericLabels=TRUE,maxBlockSize=5000,randomSeed=seed,nThreads=4,
                         saveTOMs=FALSE,verbose=0)
 pool <- intersect(deg,top[net$colors!=0])
 train <- xt[va$sample_barcode,pool,drop=FALSE]
 yy <- Surv(va$OS_time,va$OS_event)
 logmsg("Cox screening",r,k,length(pool))
 p <- apply(train,2,function(g) {
   fit <- coxph(yy~g,ties="efron")
   if(!is.finite(coef(fit))) return(NA_real_)
   pchisq(as.numeric(coef(fit)^2/vcov(fit)),1,lower.tail=FALSE)
 })
 use <- names(p)[which(p<.05)]
 inner <- integer(nrow(va))
 set.seed(seed)
 for(e in 0:1){ idx <- which(va$OS_event==e);inner[idx] <- sample(rep_len(1:10,length(idx))) }
 lp <- rep(0,nrow(te)); beta <- numeric(); lambda <- NA_real_; center <- 0
 ct <- cnt[rownames(c0),te$sample_barcode,drop=FALSE]
 testsf <- sf(ct)
 if(length(use)>1){
   xx <- train[,use,drop=FALSE]
   cv <- cv.glmnet(xx,yy,family="cox",alpha=1,standardize=TRUE,foldid=inner,type.measure="deviance")
   lambda <- cv$lambda.min
   cf <- as.matrix(coef(cv,s="lambda.min"))[,1]; beta <- cf[cf!=0]
   if(length(beta)){
    center <- mean(as.numeric(xx[,names(beta),drop=FALSE]%*%beta))
    heldx <- t(trans(sweep(ct[names(beta),,drop=FALSE],2,testsf,"/")))
    lp <- as.numeric(heldx%*%beta)-center
   }
 } else if(length(use)==1) stop("One screened gene: prespecify a fallback before resuming")
 stopifnot(all(is.finite(lp)))
 list(signature=stamp,repeat_id=r,fold_id=k,seed=seed,
      discovery=tr[,.(patient_barcode,sample_barcode,sample_type)],
      inner=data.table(patient=va$patient_barcode,inner_fold=inner),
      predictions=data.table(patient_barcode=te$patient_barcode,repeat_id=r,fold_id=k,lp=lp),
      training_genes=rownames(c0),log_geomeans=lg,size_factors=sizeFactors_from_counts<-sf(c0),
      held_size_factors=testsf,dispersion=ab,
      module=setNames(net$colors,top),DEG=deg,candidate=pool,screen_p=p,lambda=lambda,
      beta=beta,beta_symbols=sym[names(beta)],center=center,session=sessionInfo())
}
for(i in seq.int(args[1],args[2])){
 r <- as.integer((i-1)%/%10+1);k <- as.integer((i-1)%%10+1)
 path <- sprintf("checkpoints/r%02df%02d.rds",r,k)
 if(file.exists(path)){stopifnot(identical(readRDS(path)$signature,stamp));next}
 fit <- fitfold(r,k)
 tmp <- paste0(path,".tmp")
 saveRDS(fit,tmp);stopifnot(file.rename(tmp,path))
 logmsg("complete",r,k,"n_test",nrow(fit$predictions),"n_beta",length(fit$beta))
}
