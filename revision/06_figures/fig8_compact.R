src <- Sys.getenv("OSCC_DATA")
readtab <- function(step, name) read.csv(file.path(src, paste0("GSE103322_R_Step",step), name))
paired <- readtab("04","01_main_nonMaxima_paired_prefixes.csv")
tests <- readtab("04","04_focused_paired_tests_PROVISIONAL.csv")
dat <- readtab("06","03_gene_by_celltype_detection.csv")
counts <- readtab("06","02_nonMaxima_prefix_celltype_counts.csv")
cells <- c("Fibroblast","Malignant","Endothelial","T cell")
cols <- c("#159A78","#CC6500","#53799D","#786B9D")
sets <- list(c("SFRP2","SFRP4","PDGFRA","FNDC1","OGN","TK1","MYBL2"),
 c("CCN1","CCL2","ZBTB10","DAB2","JUN","KYNU","CCL20"),
 c("CTSG","P4HA1","NKAIN1","LRRC38","SPSB4","SPINK6","TNFRSF4","DEFB1"))
dat$gene[dat$gene=="CYR61"] <- "CCN1"
stopifnot(nrow(paired)==8,!anyDuplicated(paired$raw_prefix),all(tests$BH_q_two_comparisons==.0078125))
for (j in seq_along(cells)) {
  s <- counts[counts$cell_type==cells[j] & counts$n_cells>0,]
  stopifnot(sum(s$n_cells)==c(524,1177,149,1031)[j],length(unique(s$raw_prefix[s$n_cells>=10]))==c(12,9,5,9)[j])
}
draw <- function(){
 layout(rbind(c(1,2,3,3),c(4,4,5,5),c(6,6,6,6)),heights=c(1,1,.17))
 par(oma=c(.05,.25,.35,.20),family="sans",fg="#29343F",col.axis="#29343F",col.lab="#29343F",
     mgp=c(1.72,.45,0),tcl=-.19,las=1)
 for (i in 1:2){
   nm <- c("Stromal5","Proliferation2")[i]
   x <- as.matrix(paired[,paste0(nm,c("_fibroblast","_malignant"))])
   par(mar=c(2.6,if(i==1)3.9 else 2.5,2.65,.65))
   plot(NA,xlim=c(.70,2.30),ylim=extendrange(x,f=.10),xaxs="i",axes=FALSE,xlab="",ylab="")
   abline(h=0,col="#D0D7DA",lty=3)
   segments(1,x[,1],2,x[,2],col="#A4AEB2",lwd=.9)
   for(j in 1:2){points(rep(j,8),x[,j],pch=16,col=cols[j],cex=.60)
       segments(j-.14,median(x[,j]),j+.14,median(x[,j]),lwd=2,col=cols[j])}
   axis(2,cex.axis=.88);axis(1,1:2,c("Fib.","Mal."),cex.axis=.95);box(bty="l")
   if(i==1) mtext("Mean gene-wise Z score",2,line=2.45,las=0,cex=.85)
   mtext(paste0(if(i==1)"A   " else "",nm),3,line=1.03,adj=0,font=2,cex=1.03)
   mtext("BH q = 0.0078",3,line=.03,adj=0,cex=.82)
 }
 for (i in 1:3){
   genes <- sets[[i]]; m <- length(genes)
   par(mar=c(2.6,4.9,2.65,.55))
   plot(NA,xlim=c(.5,4.5),ylim=c(.55,m+.45),xaxs="i",yaxs="i",axes=FALSE,xlab="",ylab="")
   a <- dat[dat$gene %in% genes & dat$cell_type %in% cells,]
   stopifnot(nrow(a)==m*4,!anyDuplicated(a[,c("gene","cell_type")]),!anyNA(a$mean_stored_expression))
   for(j in 1:4) rect(j-.46,.55,j+.46,m+.45,border=NA,col=adjustcolor(cols[j],alpha.f=.055))
   abline(h=1:m,col="#DEE3E5",lty=3)
   for(k in 1:nrow(a)){
     j <- match(a$cell_type[k],cells); y <- m+1-match(a$gene[k],genes)
     strength <- a$mean_stored_expression[k]/max(a$mean_stored_expression)
     fill <- rgb(t((1-strength)*c(246,246,246)+strength*col2rgb(cols[j])[,1]),maxColorValue=255)
     if(a$pct_detected[k]>0) symbols(j,y,circles=sqrt(a$pct_detected[k]/100)*xinch(.082),
        inches=FALSE,add=TRUE,bg=fill,fg="#59646B")
   }
   axis(2,1:m,rev(genes),tick=FALSE,cex.axis=.91);axis(1,1:4,c("Fib.","Mal.","Endo.","T cell"),tick=FALSE,cex.axis=.91)
   mtext(c("B   Stromal and proliferation genes","C   Ultra7","D   Final8 components")[i],3,line=1.03,adj=0,font=2,cex=1.03)
 }
 par(mar=rep(0,4));plot.new()
 text(.035,.52,"Detected:",adj=0,cex=.89)
 for(j in 1:3){v<-c(.1,.3,.6)[j];x<-.15+.11*(j-1);symbols(x,.52,circles=sqrt(v)*xinch(.082),inches=FALSE,add=TRUE,bg="#B0B8BC",fg="#59646B");text(x+.018,.52,paste0(v*100,"%"),adj=0,cex=.87)}
 text(.50,.52,"Shade: panel-scaled mean expression",adj=0,cex=.89)
}
pdf("Fig8_compact.pdf",width=7.2,height=5.5,pointsize=10,useDingbats=FALSE); draw();dev.off()
for(spec in list(c("Fig8_compact.png","500"),c("Fig8_compact_preview.png","150"))){
 png(spec[1],width=7.2,height=5.5,units="in",res=as.numeric(spec[2]),pointsize=10,type="cairo");draw();dev.off()
}
