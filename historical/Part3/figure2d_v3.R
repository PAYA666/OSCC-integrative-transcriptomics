project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)

net <- step5$ultra_network

ultra_order <- c(
"CCN1",
"CCL2",
"ZBTB10",
"DAB2",
"JUN",
"CCL20",
"KYNU"
)

plot_df <- data.frame(

Gene = rep(
net$Gene,
4
),

Evidence = rep(
c(
"Acquired resistance",
"Pathway convergence",
"TCGA/WGCNA",
"Network robustness"
),
each=nrow(net)
),

Score = c(

as.numeric(
net$Strict_197 &
net$Strict_111
),

pmin(
net$LeadingEdge_Evidence_Count,
3
),

as.numeric(
net$In_WGCNA_top5000
),

as.numeric(
net$Network_Interpretable
)

)

)

plot_df$Gene <- factor(
plot_df$Gene,
levels=rev(ultra_order)
)

plot_df$Evidence <- factor(
plot_df$Evidence,
levels=c(
"Acquired resistance",
"Pathway convergence",
"TCGA/WGCNA",
"Network robustness"
)
)

p <- ggplot(
plot_df,
aes(
x=Evidence,
y=Gene,
size=Score,
fill=Score
)
)+

geom_point(
shape=21,
colour="black"
)+

scale_size_continuous(
range=c(2,10)
)+

scale_fill_gradient(
low="white",
high="#B84D43"
)+

theme_classic(
base_size=16
)+

theme(
axis.text.x=
element_text(
angle=35,
hjust=1
)
)+

labs(
title=
"Multi-layer evidence supporting the frozen Ultra7 resistance program",
x=NULL,
y=NULL,
size="Evidence strength",
fill="Evidence"
)

ggsave(
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure2/Figure2D_Ultra7_dotplot.png"),
p,
width=8,
height=5,
dpi=300
)

print(p)
