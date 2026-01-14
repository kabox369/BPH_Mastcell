library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(readr)
library(scales)
library(RcppML)
library(msigdbr)
library(clusterProfiler)
library(tidyverse)
library(msigdbr)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(irGSEA)
library(AUCell)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seurat_object <- BPH_LS_mappedLabelled_EndoMT_new 

# Figure 1B
p1 <- DimPlot(
  seurat_object, 
  reduction   = "harmony.major.pcs40.umap.dims30", 
  group.by    = "primary_type", 
  label       = FALSE, 
  cols        = my_cols_major, 
  alpha       = 1, 
  pt.size     = 2.5, 
  raster.dpi  = c(2048, 2048), 
  raster      = TRUE
) + 
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    title = NULL
  )

p1 <- p1 + 
  theme(
    legend.spacing.y = unit(1.5, "cm"),  
    legend.key.size = unit(1.5, "lines"),  
    axis.title.x = element_text(size = 10, face = "bold"),  
    axis.title.y = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 18),  
    strip.text = element_text(size = 20, face = "bold"),  
    legend.position = "right",  
    legend.margin = margin(0, 20, 0, 0), 
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 5)) 
  )

ggsave("major annotation figure.pdf", p1, width = 10, height = 8, dpi = 600)

# Figure 1C
p2 <- DimPlot(
  seurat_object, 
  reduction = "harmony.major.pcs40.umap.dims30", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_minor, 
  alpha = 1, 
  pt.size = 2, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE
)+ labs(
  x = "UMAP1",
  y = "UMAP2",
  title = NULL
)

p2 <- p2 +  
  theme(
    legend.spacing.y = unit(1, "cm"),  
    legend.text = element_text(size = 14),  
    axis.title.x = element_text(size = 12, face = "bold"),  
    axis.title.y = element_text(size = 12, face = "bold"),  
    strip.text = element_text(size = 14, face = "bold"),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
  )
print(p2)
ggsave("secondary annotation figure.pdf", p2, width = 20, height = 15, dpi = 800)

#Figure 1D
group_col <- "primary_type"

cell_order <- c("Mast cell", "Myeloid cell", "NK", "B cell", "T cell",
                "Endothelia", "Fibroblast", "Epithelia")

stopifnot(group_col %in% colnames(seurat_object@meta.data))

seurat_object[[group_col]][,1] <- factor(seurat_object[[group_col]][,1], levels = cell_order)

genes_use <- c("ACPP","KLK2","KRT5","DCN","COL1A1","LUM","VWF","SELE","PECAM1",
               "CD3D","CD3E","CD3G","MS4A1","CD79A","JCHAIN","GNLY","NKG7","TYROBP",
               "CD14","FCGR3A","LYZ","CPA3","TPSAB1","KIT")
genes_use <- intersect(genes_use, rownames(seurat_object))

p <- DotPlot(
  object    = seurat_object,
  features  = genes_use,
  group.by  = group_col,
  cols      = c("#FDE0C5", "#D7301F"),
  dot.scale = 6
) +
  scale_y_discrete(limits = rev(cell_order), drop = FALSE) +
  labs(x = NULL, y = NULL,
       color = "Average\nExpression",
       size  = "Percent\nExpressed") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line   = element_line(color = "black")
  )

p
ggsave(
  "Dot plot Major Annotation.pdf", 
  plot = p
)

#Figure 1E
seurat_object <- BPH_SLT_mappedLabelled
meta <- seurat_object@meta.data
sample_list <- unique(meta$samplename)
celltype_list <- unique(meta$primary_type)
group_list <- unique(meta$group)
cellcount <- cell_count(meta)

write.csv(cellcount, file = "cell_counts_primary.csv")

proportion <- apply(cellcount, 2, function(x) x / sum(x))
proportion <- as.data.frame(t(proportion))


proportion <- proportion %>% 
  dplyr::mutate(Group = str_split_fixed(rownames(proportion), '_', n=2)[,1]) %>% 
  dplyr::mutate(Name = str_split_fixed(rownames(proportion), '_', n=2)[,2])


proportion_plot <- gather(proportion, Celltype, Value, -Name, -Group)
proportion_plot$Value <- round(proportion_plot$Value * 100, 2)

proportion_plot$Group <- factor(proportion_plot$Group, levels = c("Large", "Small"))

write.csv(proportion, file = "cell_proportion_primary.csv")

unique(proportion_plot$Celltype)

group_order <- c("Large", "Small")  

desired_primary_order <- c("T cell","NK","B cell","Myeloid cell","Mast cell","Epithelia","Endothelia","Fibroblast")

legend_order <- desired_primary_order  


proportion_plot$Group <- factor(proportion_plot$Group,levels = group_order)

proportion_plot$Celltype <- factor(proportion_plot$Celltype,levels = legend_order)

my_cols_major <- my_cols_major[legend_order]

group_proportion <- proportion_plot %>%
  group_by(Group, Celltype) %>%
  summarise(Value = sum(Value), .groups = "drop") %>%  
  group_by(Group) %>%
  mutate(Total = sum(Value),
         Percent = Value/Total*100) %>%  
  ungroup()


p <- ggplot(group_proportion, 
            aes(x = Group, y = Percent, fill = Celltype)) +
  geom_col(position = "stack", width = 0.7, 
           color = "white", linewidth = 0.2, key_glyph = "rect") +
  coord_flip() +
  scale_fill_manual(values = my_cols_primary, 
                    breaks = legend_order,
                    guide = guide_legend(
                      reverse = TRUE,
                      keywidth = 0.5, 
                      keyheight = 0.8   
                    )) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = scales::label_number(suffix = "%")) +
  labs(
    title = "Cell Type Composition Across All Samples",
    x = NULL,
    y = NULL,
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(color = "black", size = 10, family = "Arial"),
    axis.text.y = element_text(color = "black", size = 12, face = "bold", margin = margin(r = 10)),  
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),  
    legend.spacing.y = unit(0.05, "cm"),  
    legend.title = element_text(face = "bold", size = 9),  
    legend.text = element_text(size = 8),  
    legend.margin = margin(0, 0, 0, -5),  
    plot.margin = margin(1, 1.5, 1, 1, "cm"), 
    axis.line.x = element_line(color = "black", linewidth = 0.5)
  )

ggsave(
  "Cell Type Composition Across All Samples.pdf", 
  plot = p,
  width = 18,  
  height = 6,   
  device = cairo_pdf,
  dpi = 300
)

#Figure 1G, J 
object.list <- list(Large = BPH_Large_EndoMT_Cellchat, 
                    Small = BPH_Small_EndoMT_Cellchat )
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),color.use = c("#2CA02C","#FF7F0E"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",,color.use = c("#2CA02C","#FF7F0E"))
gg1 + gg2

rankNet(cellchat, mode = "comparison", measure = "weight", stacked = T, signaling = c("ICAM","VCAM","VEGF","TGFb"),do.flip =T,cutoff.pvalue = 0,
        comparison = c(1, 2),color.use = c("#2CA02C","#FF7F0E"))

#Figure 1I
object_fibro <- subset(seurat_object,subset = primary_type = "Fibroblast")

object_Fibro_gsea_H <- irGSEA.score(object = object_fibro, assay = "RNA", 
                                    slot = "data", seeds = 123, ncores = 8,
                                    min.cells = 1, min.feature = 0,
                                    custom = F, geneset = NULL, msigdb = T, 
                                    species = "Homo sapiens", category = "H",  
                                    subcategory = NULL, geneid = "symbol",
                                    method = c("singscore"),
                                    aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                                    kcdf = 'Gaussian')
result.dge_object_Fibro_gsea_H <- irGSEA.integrate(object = object_Fibro_gsea_H, 
                                                   group.by = "group",
                                                   metadata = NULL, col.name = NULL,
                                                   method = c("singscore"))
Fibro_gsea_GO_irGSEA.heatmap.plot <- irGSEA.heatmap(object =  result.dge_object_Fibro_gsea_H, 
                                                    method = c("singscore"),
                                                    top = 20,
                                                    heatmap.width = 14,
                                                    heatmap.heigh = 4.8
)

pdf("Fibro_gsea.heatmap.plotGO.pdf", width = 12,height = 8)
print(Fibro_gsea_GO_irGSEA.heatmap.plot )
dev.off()

#Figure 1 M-N
Idents(BPH_LS_mappedLabelled) <- "primary_type"
unique(BPH_LS_mappedLabelled$primary_type)
object_Fibroblast <- subset(BPH_LS_mappedLabelled, idents = c("Fibroblast"))
object_Fibroblast  <- JoinLayers(object_Fibroblast)
object_Fibroblast [["RNA"]] <- split(object_Fibroblast [["RNA"]], f = object_Fibroblast$samplename)
object_Fibroblast <- NormalizeData(object_Fibroblast)
object_Fibroblast <- FindVariableFeatures(object_Fibroblast)
object_Fibroblast <- ScaleData(object_Fibroblast)
object_Fibroblast <- RunPCA(object_Fibroblast)
object_Fibroblast <- IntegrateLayers(object = object_Fibroblast, method = HarmonyIntegration,
                                     orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                     verbose = FALSE
)
object_Fibroblast <- RunUMAP(object_Fibroblast, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_Fibroblast <- FindNeighbors(object_Fibroblast, reduction = "harmony.pcs30", dims = 1:30)
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_Fibroblast <- RunUMAP(object_Fibroblast, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_Fibroblast <- FindNeighbors(object_Fibroblast, reduction = "harmony.pcs30", dims = 1:20)
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_Fibroblast <- IntegrateLayers(
  object = object_Fibroblast, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_Fibroblast <- RunUMAP(object_Fibroblast, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_Fibroblast <- FindNeighbors(object_Fibroblast, reduction = "harmony.pcs40", dims = 1:40)
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_Fibroblast <- RunUMAP(object_Fibroblast, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_Fibroblast <- FindNeighbors(object_Fibroblast, reduction = "harmony.pcs40", dims = 1:30)
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_Fibroblast <- FindClusters(object_Fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

object_Fibroblast@reductions$UMAP <- CreateDimReducObject(
  embeddings = object_Fibroblast@reductions$harmony.pcs40.umap.dims40@cell.embeddings[,1:2],
  key = "UMAP_", 
  assay = "RNA"
)
saveRDS(object_Fibroblast,"object_Fibroblast.rds",compress = F)
my_cols_group <- c("Small" = "#e6d153","Large" = "#9e6fab")

p1 <- DimPlot(
  object_Fibroblast, 
  reduction = "UMAP", 
  group.by = "group", 
  label = FALSE, 
  cols = my_cols_group, 
  alpha = 0.6, 
  pt.size = 0.2, 
  raster = F
)+ 
  ggtitle(NULL)
ggsave("Fibroblast UMAP figure by group.pdf", p1, width = 8, height = 6,dpi = 800)
p2 <- DimPlot(
  object_Fibroblast , 
  reduction = "UMAP", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_minor, 
  alpha = 0.6, 
  pt.size = 0.2, 
  raster = F
)+ 
  ggtitle ( NULL )
ggsave("Fibroblast cell secondary annotation figure.pdf", p2, width = 10, height = 8)

#Figure O 
unique(object_fibro@meta.data$secondary_type)
Idents(object_fibro) <- "secondary_type"

# Define gene list
custom_order <- c(
  "Fibroblast_SFRP2",
  "Myofibroblast_RGS5",
  "Stromal cell_PECAM1+ACTA2+"
)
gene_list <- c("DCN","COL1A1","FBLN1",
               "ACTA2","RGS5",
               "PECAM1","VWF","CXCR4"
)

plot_df <- FetchData(object_fibro, vars = c(gene_list, "secondary_type"), layer = "data") %>%
  mutate(
    secondary_type = factor(secondary_type, levels = custom_order)  
  ) %>%
  group_by(secondary_type) %>%
  summarise_all(function(x) mean(x = expm1(x = x))) %>%
  filter(!is.na(secondary_type)) %>% 
  arrange(secondary_type) %>%  
  column_to_rownames("secondary_type") %>%
  log2()
col_fun <- circlize::colorRamp2(
  c(-4, 0, 4),
  c("#0F7B9F", "white", "#D83215")
)
plot_df_t <- plot_df %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>% 
  mutate(Gene = factor(Gene, levels = gene_list)) %>%
  column_to_rownames("Gene") %>% 
  as.matrix()

hm <- ComplexHeatmap::Heatmap(
  plot_df_t, 
  col = col_fun,
  name = "Gene expression",
  
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 12, fontface = "italic", col = "black"),
  column_names_gp = gpar(fontsize = 12),
  
  rect_gp = gpar(col = "gray90", lwd = 0.5), 
  border_gp = gpar(col = "black", lwd = 1.5),  
  
  heatmap_legend_param = list(
    legend_direction = "vertical",
    title_position = "topcenter",
    legend_height = unit(2, "cm"),
    legend_width = unit(1, "cm"),
    gap = unit(0.5, "cm"),
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  right_annotation = rowAnnotation(
    spacer = anno_empty(
      width = unit(2, "cm"),  
      border = FALSE
    )
  ),
  width = ncol(plot_df_t)*unit(12, "mm"),  
  height = nrow(plot_df_t)*unit(10, "mm")  
)

pdf("fibro_vertical_heatmap.pdf", 
    width = 10,   
    height = 16,
    useDingbats = FALSE)
draw(hm,
     padding = unit(c(2, 6, 2, 2), "cm"),  
     heatmap_legend_side = "right",
     gap = unit(3, "cm"),  
     adjust_annotation_extension = TRUE)

dev.off()