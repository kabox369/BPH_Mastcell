#Figure S8A
suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(Matrix)
  library(scales)
  library(patchwork)
})

# === same plotting helper as before ===
plot_gene_umap_log10 <- function(cds, gene, assay_name = "counts",
                                 lim = c(0, 2),
                                 pt_size = 0.6) {
  um <- as.data.frame(reducedDims(cds)$UMAP)
  colnames(um) <- c("UMAP_1", "UMAP_2")
  um$cell <- rownames(um)
  
  gshort <- rowData(cds)$gene_short_name
  rn <- rownames(cds)
  
  if (!is.null(gshort) && gene %in% gshort) {
    ridx <- which(gshort == gene)[1]
    rnm <- rn[ridx]
  } else if (gene %in% rn) {
    rnm <- gene
  } else {
    stop("Gene not found in cds: ", gene)
  }
  
  if (!assay_name %in% names(assays(cds))) {
    stop("Assay not found in cds assays: ", assay_name,
         ". Available: ", paste(names(assays(cds)), collapse = ", "))
  }
  
  expr_mat <- assays(cds)[[assay_name]]
  v <- as.numeric(expr_mat[rnm, colnames(cds)])
  df <- data.frame(cell = colnames(cds), expr = v, stringsAsFactors = FALSE)
  df <- merge(df, um, by = "cell", all.x = TRUE)
  
  df$expr_log10 <- log10(df$expr + 1)
  df$expr_log10[df$expr == 0] <- NA
  
  ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "grey80", size = pt_size) +
    geom_point(
      data = df[!is.na(df$expr_log10), ],
      aes(color = expr_log10),
      size = pt_size
    ) +
    scale_color_viridis_c(
      name = "log10(Expression)",
      limits = lim,
      oob = squish
    ) +
    labs(title = gene, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5),
      legend.position = "right"
    )
}

# === genes in your panel ===
genes_panel <- c("PECAM1","ACTA2","ENG",
                 "MCAM","CD34","KDR",
                 "CDH5","PTPRC","CD14")

# === optional per-gene color limits (to mimic your example) ===
# You can comment this out and use a single global limit if you prefer.
lims_map <- list(
  PECAM1 = c(0, 1.2),
  ACTA2  = c(0, 2.0),
  ENG    = c(0, 2.0),
  MCAM   = c(0, 2.0),
  CD34   = c(0, 2.0),
  KDR    = c(0, 2.0),
  CDH5   = c(0, 2.0),
  PTPRC  = c(0, 2.0),
  CD14   = c(0, 2.0)
)

# build plots
plots <- lapply(genes_panel, function(g) {
  lim <- lims_map[[g]]
  plot_gene_umap_log10(cds, g, assay_name = "counts", lim = lim, pt_size = 0.55) +
    theme(legend.position = "none")
})

# collect legends: only keep one (right side), like your figure
p_last <- plot_gene_umap_log10(cds, genes_panel[1], assay_name = "counts",
                               lim = lims_map[[genes_panel[1]]], pt_size = 0.55)

# layout 3 x 3
p_grid <- wrap_plots(plots, ncol = 3) &
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text()
  )

# add a single legend on the right
# (use the first plot's legend; all plots share same scale name)
p_final <- p_grid + plot_layout(guides = "collect") & theme(legend.position = "right")

p_final

#Figure S8B-C
#need run FigS7code
DotPlot(BPH_GEO_EndMT,features = c('VWF','PECAM1','ACTA2',
                                   'MCAM','CD34','KDR',
                                   'PTPRC','CD14'), min.cutoff = 0.5, cols = c("gray", "coral2"))
FeaturePlot(BPH_GEO_EndMT,features = c('VWF','PECAM1','ACTA2',
                                       'MCAM','CD34','KDR',
                                       'PTPRC','CD14'), cols = c("gray", "coral2"),
            min.cutoff = 0.5, reduction = "umap.harmony",pt.size = 1)