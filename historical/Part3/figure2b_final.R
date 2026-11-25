project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)
library(dplyr)

plot_df <- leading_summary

plot_df$Pathway_short <- c(

"TNFA/NFKB",

"Hypoxia",

"EMT",

"TGF-beta"

)

plot_df$Pathway_short <- factor(
  plot_df$Pathway_short,
  levels=c(
    "TGF-beta",
    "EMT",
    "Hypoxia",
    "TNFA/NFKB"
  )
)

p <- ggplot(

plot_df,

aes(

x=LeadingEdge_Intersection,

y=Pathway_short

)

)+

geom_col(

fill="#B84D43",

width=0.65

)+

geom_text(

aes(
label=LeadingEdge_Intersection
),

hjust=-0.25,

size=4,

fontface="bold"

)+

scale_x_continuous(

expand=
expansion(
mult=c(0,0.15)
)

)+

labs(

title=
"Shared leading-edge genes across resistance pathways",

subtitle=
"Genes recovered in both acquired cisplatin-resistance models",

x=
"Leading-edge overlap genes",

y=NULL

)+

theme_classic(

base_size=12

)+

theme(

plot.title=
element_text(
face="bold",
size=13
),

plot.subtitle=
element_text(
size=9,
colour="grey40"
),

axis.text.y=
element_text(
size=10
)

)

print(p)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure2")

ggsave(

file.path(
out_dir,
"Figure2B_Shared_leading_edge_overlap.png"
),

p,

width=95,

height=65,

units="mm",

dpi=600,

bg="white"

)

cat("\n========== FIGURE2B RETURN ==========\n")

print(plot_df)
