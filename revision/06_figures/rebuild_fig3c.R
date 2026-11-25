library(ggplot2)

src <- paste0(Sys.getenv("OSCC_DATA"), "/Part3_Cisplatin/Step3_Gene_Convergence/Part3_Gene_Evidence_Master.tsv")
d <- read.delim(src)
x <- subset(d, Strict_197 & Strict_111)
stopifnot(nrow(x) == 129L, !anyDuplicated(x$Gene),
          all(is.finite(x$logFC_197) & is.finite(x$logFC_111)),
          all(x$logFC_197 != 0 & x$logFC_111 != 0))
group <- ifelse(x$logFC_197 * x$logFC_111 < 0, "Opposite",
                ifelse(x$logFC_197 > 0, "Up in both", "Down in both"))
counts <- as.data.frame(table(factor(group, levels=c("Down in both","Opposite","Up in both"))))
names(counts) <- c("direction","genes")
stopifnot(identical(counts$genes, c(43L,32L,54L)))
write.table(counts, "Fig3C_counts.tsv", sep="\t", row.names=FALSE, quote=FALSE)
p <- ggplot(counts, aes(direction, genes)) +
  geom_col(width=.65, fill="#B34E43") +
  geom_text(aes(label=genes), vjust=-.38, family="Arial", size=3.45) +
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,20), expand=c(0,0)) +
  labs(x=NULL,y="Number of genes",title="Conserved DEG direction") +
  theme_classic(base_family="Arial",base_size=8.2) +
  theme(plot.title=element_text(face="bold",size=11,margin=margin(b=5)),
        axis.text=element_text(colour="black"),axis.text.x=element_text(angle=30,hjust=1),
        axis.line=element_line(linewidth=.35),plot.margin=margin(3,5,1,2))
for (k in c(1,5502/3100)) {
  png(if(k==1) "Fig3C_word.png" else "Fig3C_print.png",
      width=round(875*k),height=round(744*k),res=round(260*k),type="cairo")
  print(p)
  dev.off()
}
