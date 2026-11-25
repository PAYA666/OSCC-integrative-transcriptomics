project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(reshape2)

evidence_mat <- data.frame(

Gene = c(
"CCN1",
"CCL2",
"ZBTB10",
"DAB2",
"JUN",
"KYNU",
"CCL20"
),

DEG = c(3,3,3,3,3,3,3),

Meta = c(3,3,3,3,3,3,3),

Pathway = c(3,2,2,2,3,2,2),

WGCNA = c(3,3,3,3,1,1,1),

TCGA = c(2,2,2,2,2,1,1),

GDSC1 = c(1,1,0,1,1,1,0)

)

rownames(evidence_mat) <- evidence_mat$Gene

evidence_mat$Gene <- NULL

plot_df <- melt(
  as.matrix(evidence_mat)
)

colnames(plot_df) <- c(
  "Gene",
  "Evidence",
  "Score"
)

plot_df$Gene <- factor(
  plot_df$Gene,
  levels = rev(rownames(evidence_mat))
)

p <- ggplot(
  plot_df,
  aes(
    Evidence,
    Gene,
    fill = Score
  )
)+
geom_tile(
  color="white",
  linewidth=0.5
)+
scale_fill_gradient(
  low="#F7F7F7",
  high="#9E2F25",
  limits=c(0,3)
)+
labs(
 title="Multi-layer evidence supporting frozen Ultra7 resistance candidates",
 subtitle="Integrated evidence from discovery, pathway, tumor-context and pharmacogenomic layers",
 x=NULL,
 y=NULL,
 fill="Evidence\nstrength"
)+
theme_classic(base_size=16)+
theme(
 plot.title=element_text(
  face="bold",
  size=22
 ),
 axis.text.x=element_text(
  angle=45,
  hjust=1
 )
)

ggsave(
 "Figure5A_Ultra7_multilayer_evidence.pdf",
 p,
 width=8,
 height=6
)

p
