#Figure 3A
pathways.show2 <- c("TGFb","VEGF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(BPH_EndoMT_Cellchat, signaling = pathways.show2, layout = "circle", signaling.name = c("Weighted Integration of TGF-β and VEGF Signaling Pathways"),
                    targets.use = "Stromal cell_PECAM1+ACTA2+",vertex.weight.max = 1)

object.list <- list(Large = BPH_Large_EndoMT_Cellchat, 
                    Small = BPH_Small_EndoMT_Cellchat )
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Mast cell', targets.use = 'Stromal cell_PECAM1+ACTA2+', stacked =T, do.stat = TRUE,
        signaling = c("VEGF","TGFb"),x.rotation = 0, do.flip = F,
        comparison = c(1, 2),color.use = c("#2CA02C","#FF7F0E"),)

#Figure 3B
library(ggplot2)
library(Seurat)
library(dplyr)
object_mast <- subset(BPH_LS_mappedLabelled_noEndoMT)
object_mast <- NormalizeData(object_mast)
object_mast <- FindVariableFeatures(object_mast)
object_mast <- ScaleData(object_mast)
object_mast <- RunPCA(object_mast)
object_mast <- IntegrateLayers(object = object_mast, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                               verbose = FALSE
)
object_mast <- RunUMAP(object_mast, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_mast <- FindNeighbors(object_mast, reduction = "harmony.pcs30", dims = 1:30)
object_mast <- FindClusters(object_mast, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_mast <- FindClusters(object_mast, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_mast <- FindClusters(object_mast, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_mast <- FindClusters(object_mast, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_mast <- RunUMAP(object_mast, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_mast <- FindNeighbors(object_mast, reduction = "harmony.pcs30", dims = 1:20)
object_mast <- FindClusters(object_mast, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_mast <- FindClusters(object_mast, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_mast <- FindClusters(object_mast, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_mast <- FindClusters(object_mast, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_mast <- IntegrateLayers(
  object = object_mast, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_mast <- RunUMAP(object_mast, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_mast <- FindNeighbors(object_mast, reduction = "harmony.pcs40", dims = 1:40)
object_mast <- FindClusters(object_mast, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_mast <- FindClusters(object_mast, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_mast <- FindClusters(object_mast, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_mast <- FindClusters(object_mast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_mast <- RunUMAP(object_mast, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_mast <- FindNeighbors(object_mast, reduction = "harmony.pcs40", dims = 1:30)
object_mast <- FindClusters(object_mast, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_mast <- FindClusters(object_mast, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_mast <- FindClusters(object_mast, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_mast <- FindClusters(object_mast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")


saveRDS(object_mast,file = "object_mast.rds", compress = FALSE)


marker_list_mast <- c(
  "CD3D", "CD3E", "CD3G",          # T
  "CD8A", "CD8B",                  # CD8T
  "CD40LG", "CD4",                # CD4T
  "CD79A","MZB1",                 # B
  "NKG7", "KLRF1",                # NK
  "CD14", "FCGR3A",                # macro
  "LILRA4", "IL3RA",              # pDC
  "CD1C", "CD1E", "FCER1A",        # mDC                        
  "FCER1A", "HDC",                # basophil
  "LTF","CSF3R",                  # neutrophils
  "EPCAM", "ACPP",                # epi
  "VWF", "PECAM1",                # endo
  "ACTA2", "MYH11",                # SMC
  "DCN", "LUM",                   # fibro
  "KIT","CPA3","MS4A1", "TPSAB1","TPSB2","CDK15", "GATA2"#mast
)

object_mast <- BPH_SLN.secondary.annotation.mastcell_reumap
mastcelljson <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/nameing json/mast cell.json")
object_mast <- label_clusters(object_mast,mastcelljson)
DimPlot(object = object_mast,reduction ="harmony.pcs30.umap.dims30",group.by = "secondary_type" ,label = T )
object_mast<- subset(object_mast, subset = secondary_type != "delete")
saveRDS(object_mast,file = "object_mast.rds", compress = FALSE)
object_mast@reductions$UMAP <- CreateDimReducObject(
  embeddings = object_mast@reductions$harmony.pcs30.umap.dims30@cell.embeddings[,1:2],
  key = "UMAP_", 
  assay = "RNA"
)

p1 <- DimPlot(
  object_mast , 
  reduction = "UMAP", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_secondary, 
  alpha = 0.5, 
  pt.size = 1, 
  raster = F
)

# 调整图形样式
p1 <- p1  +  
  theme(
    legend.spacing.y = unit(2, "cm"),  
    legend.key.size = unit(1.5, "lines"),  
    legend.key.height = unit(1.5, "cm"),  
    axis.title.x = element_text(size = 10, face = "bold"),  
    axis.title.y = element_text(size = 10, face = "bold"), 
    legend.text = element_text(size = 22),
    strip.text = element_text(size = 20, face = "bold"), 
    legend.position = "right",  
    legend.direction = "vertical",  
    legend.margin = margin(0, 20, 0, 0),  
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 6))  
  )

ggsave("Mast Cellsecondary annotation figure.pdf", p1, width = 8, height = 6)

#Figure 3C
p6 <- DimPlot(
  object_Mast , 
  reduction = "UMAP", 
  group.by = "group", 
  label = FALSE, 
  cols = my_cols_group, 
  alpha = 0.8, 
  pt.size = 1, 
  raster = F
)

p2 <- p2 + 
  theme(
    legend.spacing.y = unit(2, "cm"),  
    legend.key.size = unit(1.5, "lines"),  
    legend.key.height = unit(1.5, "cm"),  
    axis.title.x = element_text(size = 10, face = "bold"),  
    axis.title.y = element_text(size = 10, face = "bold"),  
    legend.text = element_text(size = 22), 
    strip.text = element_text(size = 20, face = "bold"),  
    legend.position = "right",  
    legend.direction = "vertical",
    legend.margin = margin(0, 20, 0, 0), 
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 6))  
  )
ggsave("Mast cell UMAP by group.pdf", p2, width = 12, height = 10)

#Figure 3D
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)  # install if needed: remotes::install_github("satijalab/seurat-wrappers")
})

reduction <- "harmony.pcs30.umap.dims30"
cds <- as.cell_data_set(object_mast)
umap_mat <- Embeddings(object_mast, reduction = reduction)
# Ensure same cell order
umap_mat <- umap_mat[colnames(cds), , drop = FALSE]
reducedDims(cds)$UMAP <- umap_mat

# If your clusters are stored elsewhere, replace "seurat_clusters" accordingly.
colData(cds)$seurat_clusters <- Idents(object_mast)[colnames(cds)] |> as.character()

cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(cds, use_partition = TRUE)

# Option A: root by a cluster id (edit "0" to the cluster you want as start)
root_cluster <- "Mast cell_JAK2"
root_cells <- colnames(cds)[colData(cds)$seurat_clusters == root_cluster]

cds <- order_cells(cds, root_cells = root_cells)

plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_groups_by_cluster = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)