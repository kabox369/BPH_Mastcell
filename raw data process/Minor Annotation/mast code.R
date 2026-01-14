unique(object_merged$primary_type)
library(rjson)
Idents(object = object_merged) <- "primary_type"
object_Mast <- subset(object_merged,idents = "Mast cell")
object_Mast <- NormalizeData(object_Mast)
object_Mast <- FindVariableFeatures(object_Mast)
object_Mast <- ScaleData(object_Mast)
object_Mast <- RunPCA(object_Mast)
object_Mast <- IntegrateLayers(object = object_Mast, method = HarmonyIntegration,
                            orig.reduction = "pca", new.reduction = "harmony.secondary.pcs30", ndim = 30,
                            verbose = FALSE
)
object_Mast <- RunUMAP(object_Mast, reduction = "harmony.secondary.pcs30", dims = 1:30, reduction.name = "harmony.secondary.pcs30.umap.dims30")
object_Mast <- FindNeighbors(object_Mast, reduction = "harmony.secondary.pcs30", dims = 1:30)
object_Mast <- FindClusters(object_Mast, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3")
object_Mast <- FindClusters(object_Mast, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6")
object_Mast <- FindClusters(object_Mast, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9")
object_Mast <- FindClusters(object_Mast, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2")

object_Mast <- RunUMAP(object_Mast, reduction = "harmony.secondary.pcs30", dims = 1:20, reduction.name = "harmony.secondary.pcs30.umap.dims20")
object_Mast <- FindNeighbors(object_Mast, reduction = "harmony.secondary.pcs30", dims = 1:20)
object_Mast <- FindClusters(object_Mast, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3")
object_Mast <- FindClusters(object_Mast, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6")
object_Mast <- FindClusters(object_Mast, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9")
object_Mast <- FindClusters(object_Mast, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2")

object_Mast <- IntegrateLayers(
  object = object_Mast, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs40", ndim = 40,
  verbose = FALSE
)
object_Mast <- RunUMAP(object_Mast, reduction = "harmony.secondary.pcs40", dims = 1:40, reduction.name = "harmony.secondary.pcs40.umap.dims40")
object_Mast <- FindNeighbors(object_Mast, reduction = "harmony.secondary.pcs40", dims = 1:40)
object_Mast <- FindClusters(object_Mast, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3")
object_Mast <- FindClusters(object_Mast, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6")
object_Mast <- FindClusters(object_Mast, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9")
object_Mast <- FindClusters(object_Mast, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2")

object_Mast <- RunUMAP(object_Mast, reduction = "harmony.secondary.pcs40", dims = 1:30, reduction.name = "harmony.secondary.pcs40.umap.dims30")
object_Mast <- FindNeighbors(object_Mast, reduction = "harmony.secondary.pcs40", dims = 1:30)
object_Mast <- FindClusters(object_Mast, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3")
object_Mast <- FindClusters(object_Mast, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6")
object_Mast <- FindClusters(object_Mast, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9")
object_Mast <- FindClusters(object_Mast, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2")

DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()

marker_list_mast <- c(
  "CD3D", "CD3E", "CD3G",          # T
  "CD8A", "CD8B",                  # CD8T
  "CD40LG", "CD4",                # CD4T
  "CD79A","MZB1",                 # B
  "NKG7", "KLRF1",                # NK
  "CD14", "FCGR3A",                # macro
  "IL3RA",              # pDC
  "CD1C", "CD1E", "FCER1A",        # mDC                        
  "FCER1A", "HDC",                # basophil
  "LTF","CSF3R",                  # neutrophils
  "EPCAM", "ACPP",                # epi
  "VWF", "PECAM1",                # endo
  "ACTA2", "MYH11",                # SMC
  "DCN", "LUM",                   # fibro
  "KIT","CPA3","MS4A1", "TPSAB1","TPSB2","CDK15", "GATA2",#mast
  "JAK2","CDKN1A","CTSG"
)



VlnPlot(object_Mast,features = marker_list_mast,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2") + NoLegend()
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "group",raster=F)
DimPlot(object_Mast, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()

Mastjson <- rjson::fromJSON(file = "~/Single Cell/BPH_LS/nameing json/Mast cell.json")
object_Mast <- label_clusters(object_Mast,Mastjson)
DimPlot(object = object_Mast,reduction ="harmony.secondary.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
VlnPlot(object_Mast,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
object_Mast<- subset(object_Mast, subset = secondary_type != "delete")
saveRDS(object_Mast,file = "BPH_LS.secondary.annotation.Mast.rds", compress = FALSE)

