library(jsonlite)
Idents(object = object_merged) <- "primary_type"
object_B <- subset(object_merged,idents = "B cell")
object_B <- NormalizeData(object_B)
object_B <- FindVariableFeatures(object_B)
object_B <- ScaleData(object_B)
object_B <- RunPCA(object_B)
object_B <- IntegrateLayers(object = object_B, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony.secondary.pcs30", ndim = 30,
                               verbose = FALSE
)
object_B <- RunUMAP(object_B, reduction = "harmony.secondary.pcs30", dims = 1:30, reduction.name = "harmony.secondary.pcs30.umap.dims30")
object_B <- FindNeighbors(object_B, reduction = "harmony.secondary.pcs30", dims = 1:30)
object_B <- FindClusters(object_B, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3")
object_B <- FindClusters(object_B, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6")
object_B <- FindClusters(object_B, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9")
object_B <- FindClusters(object_B, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2")

object_B <- RunUMAP(object_B, reduction = "harmony.secondary.pcs30", dims = 1:20, reduction.name = "harmony.secondary.pcs30.umap.dims20")
object_B <- FindNeighbors(object_B, reduction = "harmony.secondary.pcs30", dims = 1:20)
object_B <- FindClusters(object_B, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3")
object_B <- FindClusters(object_B, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6")
object_B <- FindClusters(object_B, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9")
object_B <- FindClusters(object_B, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2")

object_B <- IntegrateLayers(
  object = object_B, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs40", ndim = 40,
  verbose = FALSE
)
object_B <- RunUMAP(object_B, reduction = "harmony.secondary.pcs40", dims = 1:40, reduction.name = "harmony.secondary.pcs40.umap.dims40")
object_B <- FindNeighbors(object_B, reduction = "harmony.secondary.pcs40", dims = 1:40)
object_B <- FindClusters(object_B, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3")
object_B <- FindClusters(object_B, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6")
object_B <- FindClusters(object_B, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9")
object_B <- FindClusters(object_B, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2")

object_B <- RunUMAP(object_B, reduction = "harmony.secondary.pcs40", dims = 1:30, reduction.name = "harmony.secondary.pcs40.umap.dims30")
object_B <- FindNeighbors(object_B, reduction = "harmony.secondary.pcs40", dims = 1:30)
object_B <- FindClusters(object_B, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3")
object_B <- FindClusters(object_B, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6")
object_B <- FindClusters(object_B, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9")
object_B <- FindClusters(object_B, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2")

DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_B,file = "BPH_NSL.subclass.T_NK.rds", compress = FALSE)




marker_list_Bcell <- c(
  "CD3D", "CD3E", "CD3G",          # T
  "CD8A", "CD8B",                  # CD8T
  "CD40LG", "CD4",                # CD4T
  "NKG7", "KLRF1",                # NK
  "CD14", "FCGR3A",                # macro
  "KIT", "CPA3",                  # mast
  "LILRA4", "IL3RA",              # pDC
  "CD1C", "CD1E", "FCER1A",        # mDC                        
  "FCER1A", "HDC",                # basophil
  "LTF","CSF3R",                  # neutrophils
  "EPCAM", "ACPP",                # epi
  "VWF", "PECAM1",                # endo
  "ACTA2", "MYH11",                # SMC
  "DCN", "LUM",                    # fibro
  "IGHD", "TCL1A", "FCER2", "IGHM", # Bn
  "AIM2", "FCRL5", "CD27", "TNFRSF13B", "TXNP", "GPR183", # Bm
  "BCL6", "AICDA", "RGS13", "IL21R", # Bgc
  "JCHAIN","MZB1", "CD38", "PRDM1", "XBP1", # plasma
  "MKI67", "STMM1", "TOP2A", "HMGB2" # cycling B
)


VlnPlot(object_B,features = marker_list_Bcell,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2") + NoLegend()
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename",ncol = 3,raster=F)
DimPlot(object_B, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()

Bjson <- rjson::fromJSON(file = "~/Single Cell/BPH_LS/nameing json/B cell.json")
object_B <- label_clusters(object_B,Bjson)
DimPlot(object = object_B,reduction ="harmony.secondary.pcs40.umap.dims30",group.by = "secondary_type" ,label = T )
VlnPlot(object_B,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
object_B<- subset(object_B, subset = secondary_type != "delete")
saveRDS(object_B,file = "BPH_LS.secondary.annotation.B.rds", compress = FALSE)
\