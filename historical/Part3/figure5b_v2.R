project_root <- normalizePath(Sys.getenv("OSCC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)

library(ggplot2)

out_dir <-
file.path(project_root, "Part3_Cisplatin/Figures_FINAL/Figure5")

df <- data.frame(

y = 5:1,

label = c(
"Acquired cisplatin-resistance models\nTwo independent bulk datasets",
"Cross-dataset convergence\ngene and pathway architecture",
"Human tumor projection\nTCGA and WGCNA context",
"Orthogonal pharmacogenomic assessment\nGDSC1 HNSCC",
"Frozen Ultra7 resistance program"
)

)

p <- ggplot(
df,
aes(
x=0,
y=y
)
)+

geom_label(
aes(
label=label
),

size=5,

label.size=0.6,

fill="white",

fontface="bold",

width=0.55

)+

geom_segment(
data=data.frame(
x=0,
xend=0,
y=4.65:1.65,
yend=4.35:1.35
),

aes(
x=x,
xend=xend,
y=y,
yend=yend
),

arrow=
arrow(
length=unit(0.25,"cm")
),

linewidth=0.8
)+

annotate(
"point",
x=0,
y=1,
size=10
)+

theme_void()+

labs(
title=
"Evidence hierarchy supporting the frozen Ultra7 program"
)+

theme(
plot.title=
element_text(
size=20,
face="bold",
hjust=0.5
)
)

ggsave(
file.path(
out_dir,
"Figure5B_Evidence_hierarchy_V2.pdf"
),

p,

width=7,

height=9
)

print(p)
