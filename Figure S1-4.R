library(ggplot2)
library(Seurat)
library(patchwork)
library(extrafont)
library(dplyr)
seurat_object <- BPH_LS_mappedLabelled
seurat_object@meta.data$nCount_RNA <- as.numeric(seurat_object@meta.data$nCount_RNA)
seurat_object@meta.data$nFeature_RNA <- as.numeric(seurat_object@meta.data$nFeature_RNA)
seurat_object@meta.data$percent.mt <- as.numeric(seurat_object@meta.data$percent.mt)
unique(seurat_object@meta.data$group)
colnames(seurat_object@reductions[["UMAP"]]@cell.embeddings) =c('UMAP_1','UMAP_2')
head(seurat_object@reductions[["UMAP"]]@cell.embeddings)
colnames(seurat_object@reductions[["UMAP"]]@cell.embeddings) <- c('UMAP_1', 'UMAP_2')

if (is.null(seurat_object@reductions[["UMAP"]])) {
  seurat_object <- RunUMAP(seurat_object,reduction = "pca", dims = 1:30 ) 
}
seurat_object@reductions$harmony.major.pcs40.umap.dims30

#Figure S1
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
cols     <- c("#ED9F9B", '#d4de9c', '#b7deea')

plots <- lapply(seq_along(features), function(i) {
  VlnPlot(
    seurat_object,
    features = features[i],
    group.by = "samplename",
    flip     = TRUE,
    pt.size  = 0,
    raster   = TRUE
  ) +
    scale_fill_manual(values = rep(cols[i], length(unique(seurat_object$samplename)))) +
    theme_minimal(base_size = 16) +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold", size = 15),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y  = element_text(size = 12),
      panel.grid   = element_blank(),
      axis.title   = element_blank(),
      plot.margin  = margin(10, 10, 10, 10)
    ) +
    NoLegend() +
    ggtitle(features[i])
})

combined <- wrap_plots(plots, ncol = 1)
ggsave("QC_pic_sample.pdf", combined, width = 10, height = 8)

#Figure S2
p2 <- DimPlot(
  seurat_object, 
  reduction = "harmony.major.pcs40.umap.dims30", 
  group.by = "primary_type", 
  label = FALSE, 
  cols = my_cols_major, 
  alpha = 1, 
  pt.size = 2.5, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE,
  split.by = "samplename", 
  ncol = 5
)+ labs(
  x = "UMAP1",
  y = "UMAP2",
  title = NULL
)

p2 <- p2 + 
  theme(
    legend.spacing.y = unit(1, "cm"), 
    legend.text = element_text(size = 12), 
    axis.title.x = element_text(size = 12, face = "bold"),  
    axis.title.y = element_text(size = 12, face = "bold"),  
    strip.text = element_text(size = 14, face = "bold"), 
    panel.spacing = unit(1, "lines") 
  ) +
  facet_wrap(~samplename, ncol = 5, scales = "fixed")  

ggsave("major annotation figure by samplename.pdf", p2, width = 30, height = 15, dpi = 600)


#Figure S3
p3 <- DimPlot(
  seurat_object, 
  reduction = "harmony.major.pcs40.umap.dims30", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_secondary, 
  alpha = 1, 
  pt.size = 2.5, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE,
  split.by = "samplename", 
  ncol = 5
)+ labs(
  x = "UMAP1",
  y = "UMAP2",
  title = NULL
)

p3 <- p3  +  
  theme(
    legend.spacing.y = unit(1, "cm"),  
    legend.text = element_text(size = 12), 
    axis.title.x = element_text(size = 12, face = "bold"),  
    axis.title.y = element_text(size = 12, face = "bold"),  
    strip.text = element_text(size = 14, face = "bold"),  
    panel.spacing = unit(1, "lines")
  ) +
  facet_wrap(~samplename, ncol = 5, scales = "fixed") 
ggsave("secondary annotation figure by samplename.pdf", p3, width = 30, height = 15, dpi = 600)

#Figure S4
unique(seurat_object@meta.data$secondary_type)
"secondary_type" <- Idents(seurat_object) 
desired_secondary_order <- c("CD4+Tnaive_LEF1","CD4+Tem_NR4A2","CD4+Trm_CLNK","CD4+Treg_FOXP3","CD4+Tex_CTLA4", 
                             "CD8+Teff_GZMK","CD8+Temra_GNLY", "CD8+Trm_ZNF683","Proliferative T","CD8+MAIT_SLC4A10", 
                             "CD56+CD16- NK","CD56+CD16+ NK",
                             "B naive", "B memory", "Plasma",
                             "Macrophage_FCN1", "Macrophage_CX3CR1", "Macrophage_NR4A3", "mDC","Neutrophil",
                             "Mast cell_CDKN1A","Mast cell_CTSG", "Mast cell_JAK2",
                             "Endothelia", 
                             "Fibroblast_SFRP2","Myofibroblast_RGS5",
                             "Stromal cell_PECAM1+ACTA2+",
                             "Basal epi_KRT5", "Intermediate epi_KRT5",
                             "Luminal epi_EPCAM","Luminal epi_RPL","Luminal epi_DCN",
                             "Club cell_SCGB3A1" )  
seurat_object@meta.data$secondary_type <- factor(
  seurat_object@meta.data$secondary_type,
  levels = desired_secondary_order  
)
secondaryfeature <- c("CD3D","CD3E","CD3G",
                      "CCR7","SELL","LEF1","IL7R","TCF7",
                      "CD4","CD40LG","NR4A2","CLNK","FOXP3","PDCD1","CTLA4","LAG3",
                      "CD8A","CD8B","NKG7","GNLY","GZMB","GZMA","GZMK","GZMH",
                      "ITGAE","HAVCR2","ZNF683","MKI67",
                      "SLC4A10","TRAV1-2",
                      "NCAM1","FCGR3A","KLRF1","KLRD1",
                      "MS4A1","CD79A","CD79B",
                      "IGHD","TCL1A","FCER2","IGHM",
                      "AIM2","TNFRSF13B",
                      "JCHAIN","MZB1",
                      "LYZ","CD14","CD68","CD163","C1QA","FCN1","CX3CR1","NR4A3",
                      "KIT","CPA3","TPSAB1","GATA2","CDKN1A","CDKN2A","TGFB1","VEGFA","JAK2","CXCL8",
                      "VWF","PECAM1",
                      "DCN","COL1A1","COL1A2","ACTA2","RGS5","SFRP2",
                      "EPCAM","ACPP","KLK2",
                      "KRT13","KRT5","KRT14","KRT15","TP63",
                      "SCGB3A1"
)

P <- VlnPlot(seurat_object,features = secondaryfeature,cols = zzm60colors,raster = F,pt.size = 0,flip = T, stack = T,alpha = 0.1,group.by ="secondary_type")+ NoLegend()
ggsave("Vinplot secondary type.pdf",P,width = 25,height = 30)

