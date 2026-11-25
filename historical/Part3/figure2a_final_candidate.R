rm(list = ls())

project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
options(stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)

step2_file <-
file.path(project_root, "Part3_Cisplatin/Step2_Hallmark_GSEA_V2/Part3_Step2_Hallmark_V2_COMPLETE.rds")

step2 <- readRDS(
  step2_file
)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure2")

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

g197 <- as.data.frame(
  step2$GSE197561$all
)

g111 <- as.data.frame(
  step2$GSE111585$all
)

cat("\nGSE197561 columns:\n")

print(colnames(g197))

cat("\nGSE111585 columns:\n")

print(colnames(g111))

core_pathways <- c(

"HALLMARK_TNFA_SIGNALING_VIA_NFKB",

"HALLMARK_HYPOXIA",

"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",

"HALLMARK_TGF_BETA_SIGNALING"

)

cat("\n========== PATHWAY CHECK ==========\n")

for(p in core_pathways){

cat(
p,
":",
p %in% g197$pathway,
"|",
p %in% g111$pathway,
"\n"
)

}

plot_df <- bind_rows(

data.frame(

Dataset="GSE197561",

Pathway=g197$pathway[
match(
core_pathways,
g197$pathway
)
],

NES=g197$NES[
match(
core_pathways,
g197$pathway
)

]

),

data.frame(

Dataset="GSE111585",

Pathway=g111$pathway[
match(
core_pathways,
g111$pathway
)
],

NES=g111$NES[
match(
core_pathways,
g111$pathway
)

]

)

)

plot_df <- plot_df[
!is.na(plot_df$NES),
]

plot_df$Pathway <- factor(

plot_df$Pathway,

levels=rev(core_pathways),

labels=c(

"TGF-beta signaling",

"EMT",

"Hypoxia",

"TNFA/NFKB"

)

)

print(plot_df)

p <- ggplot(

plot_df,

aes(

x=NES,

y=Pathway,

fill=Dataset

)

)+

geom_col(

position=position_dodge(
width=0.7
),

width=0.62

)+

geom_vline(

xintercept=0,

linewidth=0.35,

colour="grey40"

)+

scale_fill_manual(

values=c(

"GSE197561"="#B84C43",

"GSE111585"="#3B6FA8"

)

)+

labs(

title=
"Conserved pathway activation in cisplatin resistance",

subtitle=
"Hallmark GSEA across independent acquired-resistance models",

x=
"Normalized enrichment score (NES)",

y=NULL,

fill=NULL

)+

theme_classic(

base_size=11

)+

theme(

plot.title=
element_text(
face="bold",
size=12
),

plot.subtitle=
element_text(
size=9,
colour="grey40"
),

legend.position="top",

axis.text.y=
element_text(
size=10
)

)

print(p)

ggsave(

file.path(
out_dir,
"Figure2A_Hallmark_NES_conserved_pathways.png"
),

p,

width=110,

height=75,

units="mm",

dpi=600,

bg="white"

)

ggsave(

file.path(
out_dir,
"Figure2A_Hallmark_NES_conserved_pathways.pdf"
),

p,

width=110,

height=75,

units="mm"

)

cat("\n========== FIGURE2A RETURN ==========\n")

print(plot_df)
