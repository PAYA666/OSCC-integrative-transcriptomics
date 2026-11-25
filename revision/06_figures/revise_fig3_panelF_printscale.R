library(data.table)
library(ggplot2)
library(openxlsx)
library(magick)

source_figure <- "Fig3_source.png"
source_tiff <- paste0(Sys.getenv("OSCC_DATA"), "/Integration_Manuscript/Figure3_CisplatinResistance/Figure3_CisplatinResistance_V1_600dpi.tiff")
source_table <- paste0(Sys.getenv("OSCC_DATA"), "/Integration_Manuscript/biaogewenjian/Supplementary_Tables_FINAL/Supplementary_Table_S12_Priority52_Ultra7.xlsx")

raw <- read.xlsx(source_table, sheet = "Priority52_Evidence", colNames = FALSE)
column_names <- as.character(unlist(raw[2, ], use.names = FALSE))
evidence <- raw[-c(1, 2), , drop = FALSE]
names(evidence) <- column_names
evidence <- evidence[as.logical(evidence[["In_Ultra7"]]), , drop = FALSE]
stopifnot(nrow(evidence) == 7L,
          all(as.logical(evidence[["Strict_197"]])),
          all(as.logical(evidence[["Strict_111"]])),
          all(as.logical(evidence[["Same_direction"]])))

genes <- c("CCN1", "CCL2", "ZBTB10", "DAB2", "JUN", "CCL20", "KYNU")
evidence <- evidence[match(genes, evidence[["Gene"]]), , drop = FALSE]
stopifnot(!anyNA(evidence[["Gene"]]))

columns <- c(
  "DEG\n197561", "DEG\n111585", "NF-kB\nLE",
  "Hypoxia\nLE", "EMT\nLE"
)
values <- cbind(
  as.integer(as.logical(evidence[["Strict_197"]])),
  as.integer(as.logical(evidence[["Strict_111"]])),
  as.integer(as.logical(evidence[["TNFA_NFKB_LE197"]])) +
    as.integer(as.logical(evidence[["TNFA_NFKB_LE111"]])),
  as.integer(as.logical(evidence[["Hypoxia_LE197"]])) +
    as.integer(as.logical(evidence[["Hypoxia_LE111"]])),
  as.integer(as.logical(evidence[["EMT_LE197"]])) +
    as.integer(as.logical(evidence[["EMT_LE111"]]))
)
colnames(values) <- columns
rownames(values) <- genes
long <- as.data.table(as.table(values))
setnames(long, c("Gene", "Evidence", "Value"))
long[, Gene := factor(as.character(Gene), levels = rev(genes))]
long[, Evidence := factor(as.character(Evidence), levels = columns)]
long[, binary := as.character(Evidence) %chin% columns[1:2]]
long[, label := ifelse(binary, ifelse(Value == 1, "Yes", "No"), as.character(Value))]
long[, tile_fill := fifelse(binary & Value == 1, "#3B6EA8",
                      fifelse(binary, "#F1F3F5",
                      fifelse(Value == 2, "#B94D42",
                      fifelse(Value == 1, "#F2B7AD", "#FFFFFF"))))]
long[, label_color := fifelse(tile_fill == "#B94D42" | tile_fill == "#3B6EA8",
                              "#FFFFFF", "#263238")]

panel <- ggplot(long, aes(Evidence, Gene)) +
  geom_tile(aes(fill = tile_fill), color = "#D4DADF", linewidth = 0.55,
            width = 0.94, height = 0.88) +
  geom_text(aes(label = label, color = label_color), size = 3,
            fontface = "bold", family = "Arial") +
  scale_fill_identity() + scale_color_identity() +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Ultra7 selection evidence",
       subtitle = "Pathway entries count models with leading-edge support (0-2).",
       x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(plot.title = element_text(size = 10, face = "bold", color = "#1B252B",
                                 margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.5, color = "#5B666D",
                                     margin = margin(b = 2)),
        axis.text.x = element_text(size = 6.5, face = "bold", color = "#263238",
                                   lineheight = 0.95, margin = margin(b = 6)),
        axis.text.y = element_text(size = 7.5, face = "bold", color = "#263238"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(3, 5, 2, 5))

ggsave("Fig3_panelF_printscale.png", panel, width = 4.7,
       height = 1.92, units = "in", dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave("Fig3_panelF_printscale_600.png", panel, width = 4.7,
       height = 1.92, units = "in", dpi = 600,
       bg = "white", limitsize = FALSE)
ggsave("Fig3_panelF_printscale.pdf", panel, width = 4.7,
       height = 1.92, units = "in", bg = "white",
       limitsize = FALSE, device = cairo_pdf)

base <- image_read(source_figure)
panel_image <- image_read("Fig3_panelF_printscale.png")
info <- image_info(base)
stopifnot(info$width == 3100L, info$height == 1977L,
          image_info(panel_image)$width == 1410L,
          image_info(panel_image)$height == 576L)
revised <- image_composite(base, panel_image, operator = "over",
                           offset = "+1690+1400")
image_write(revised, "Fig3_revised_printscale.png", format = "png")

tiff_base <- image_read(source_tiff)
highres_panel <- image_resize(image_read("Fig3_panelF_printscale_600.png"),
                              "2502x1021!")
tiff_revised <- image_composite(tiff_base, highres_panel, operator = "over",
                                offset = "+3000+2484")
image_write(tiff_revised, "Fig3_revised_printscale_600dpi.tiff", format = "tiff",
            compression = "lzw", density = "600x600")
write.table(as.data.frame(values), "Fig3_panelF_source.tsv", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = NA)
