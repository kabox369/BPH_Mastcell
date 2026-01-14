#Figure S9A
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

obj <- BPH_LS_mappedLabelled_EndoMT
group_col <- "primary_type"

genes_ecm <- c("COL1A2","COL4A2","COL4A1","COL1A1",
               "COL6A2","COL6A1","LAMC1","FN1","LAMA4",
               "SELE","PTN")

pal <- c("#E69F00","#56B4E9","#009E73","#F0E442",
         "#0072B2","#D55E00","#CC79A7","#999999","#8C564B")

get_data_mat <- function(object, assay = NULL, layer = "data") {
  if (is.null(assay)) assay <- DefaultAssay(object)
  
  if ("LayerData" %in% getNamespaceExports("Seurat")) {
    mat <- tryCatch(
      LayerData(object = object, assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (!is.null(mat)) return(mat)
  }
  
  GetAssayData(object = object, assay = assay, slot = layer)
}

mat <- get_data_mat(obj, layer = "data")

genes_use <- intersect(genes_ecm, rownames(mat))
missing <- setdiff(genes_ecm, genes_use)
if (length(missing) > 0) message("Missing genes: ", paste(missing, collapse = ", "))

stopifnot(length(genes_use) >= 3)
stopifnot(group_col %in% colnames(obj@meta.data))

## per-gene z-score across cells, then per-cell mean (z-mean)
mat_sub <- mat[genes_use, , drop = FALSE]

mat_z <- t(scale(t(as.matrix(mat_sub))))
ecm_zmean <- Matrix::colMeans(mat_z, na.rm = TRUE)

obj$ECMscore_zmean <- ecm_zmean

## order groups by median ECMscore_zmean (descending)
ord <- FetchData(obj, vars = c("ECMscore_zmean", group_col)) |>
  group_by(.data[[group_col]]) |>
  summarise(med = median(ECMscore_zmean, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(med)) |>
  pull(.data[[group_col]])

obj[[group_col]] <- factor(obj[[group_col]][,1], levels = ord)

## plot
p <- VlnPlot(
  obj,
  features = "ECMscore_zmean",
  group.by = group_col,
  pt.size = 0
) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = NA,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = rep(pal, length.out = length(levels(obj[[group_col]])))) +
  labs(x = NULL, y = "ECMscore (z-mean)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
    axis.ticks.length = unit(3, "pt"),
    plot.margin = margin(t = 10, r = 10, b = 25, l = 10)
  )

p

#Figure S9B-C
library(circlize)
cc = BPH_Large_EndoMT_Cellchat
cells_B <- c("Epithelia", "Fibroblast")
paths_B <- c("COLLAGEN", "LAMININ", "FN1")

cols_B <- c(
  "Epithelia" = "#1f77b4",
  "Fibroblast" = "#8c564b",
  "Stromal cell_PECAM1+ACTA2+" = "#56B4E9"
)

circos.clear()
netVisual_chord_gene(
  object = cc,
  sources.use = "Stromal cell_PECAM1+ACTA2+", 
  targets.use = cells_B,
  signaling = paths_B,
  slot.name = "net",
  color.use = cols_B,
  lab.cex = 0.45,
  small.gap = 0.6,
  big.gap = 6,
  legend.pos.x = 12,
  legend.pos.y = 1,
  title.name = paste(paths_B, collapse = " / ")
)

cells_C <- c("Myeloid cell", "Mast cell", "T cell", "B cell", "NK")
paths_C <- c("COLLAGEN", "LAMININ", "FN1", "PTN", "SELE")

cols_C <- c(
  "Stromal cell_PECAM1+ACTA2+" = "#56B4E9",
  "Myeloid cell" = "#d62728",
  "Mast cell"    = "#9467bd",
  "T cell"       = "#ff7f0e",
  "B cell"       = "#e377c2",
  "NK"           = "#bcbd22"
)

circos.clear()
netVisual_chord_gene(
  object = cc,
  sources.use = "Stromal cell_PECAM1+ACTA2+", 
  targets.use = cells_C,
  signaling = paths_C,
  slot.name = "net",
  color.use = cols_C,
  lab.cex = 0.42,
  small.gap = 0.6,
  big.gap = 6,
  legend.pos.x = 12,
  legend.pos.y = 1,
  title.name = paste(paths_C, collapse = " / ")
)

