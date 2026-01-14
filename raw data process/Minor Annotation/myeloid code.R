unique(object_merged$primary_type)
library(rjson)
Idents(object = object_merged) <- "primary_type"
object_Myeloid <- subset(object_merged,idents = "Myeloid cell")
object_Myeloid <- NormalizeData(object_Myeloid)
object_Myeloid <- FindVariableFeatures(object_Myeloid)
object_Myeloid <- ScaleData(object_Myeloid)
object_Myeloid <- RunPCA(object_Myeloid)
object_Myeloid <- IntegrateLayers(object = object_Myeloid, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony.secondary.pcs30", ndim = 30,
                               verbose = FALSE
)
object_Myeloid <- RunUMAP(object_Myeloid, reduction = "harmony.secondary.pcs30", dims = 1:30, reduction.name = "harmony.secondary.pcs30.umap.dims30")
object_Myeloid <- FindNeighbors(object_Myeloid, reduction = "harmony.secondary.pcs30", dims = 1:30)
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2")

object_Myeloid <- RunUMAP(object_Myeloid, reduction = "harmony.secondary.pcs30", dims = 1:20, reduction.name = "harmony.secondary.pcs30.umap.dims20")
object_Myeloid <- FindNeighbors(object_Myeloid, reduction = "harmony.secondary.pcs30", dims = 1:20)
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2")

object_Myeloid <- IntegrateLayers(
  object = object_Myeloid, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs40", ndim = 40,
  verbose = FALSE
)
object_Myeloid <- RunUMAP(object_Myeloid, reduction = "harmony.secondary.pcs40", dims = 1:40, reduction.name = "harmony.secondary.pcs40.umap.dims40")
object_Myeloid <- FindNeighbors(object_Myeloid, reduction = "harmony.secondary.pcs40", dims = 1:40)
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2")

object_Myeloid <- RunUMAP(object_Myeloid, reduction = "harmony.secondary.pcs40", dims = 1:30, reduction.name = "harmony.secondary.pcs40.umap.dims30")
object_Myeloid <- FindNeighbors(object_Myeloid, reduction = "harmony.secondary.pcs40", dims = 1:30)
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9")
object_Myeloid <- FindClusters(object_Myeloid, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2")

DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()

marker_list_myeloid <- c(
  "CD3D", "CD3E", "CD3G",          # T
  "CD8A", "CD8B", "GZMK",          # CD8T
  "CD40LG", "CD4",                # CD4T
  "MS4A1", "CD79A",                # B
  "NKG7", "KLRF1",                # NK                   
  "EPCAM", "ACPP",                # epi
  "VWF", "PECAM1",                # endo
  "ACTA2", "MYH11",                # SMC
  "DCN", "LUM",                    # fibro
  "FCN1", "CD14", "MX1", "S100A8", "S100A9", "S100A12", "VCAN",  # CD14+mono(经典)
  "FCGR3A", "LST1", "LILRB2", "CDKN1C", "TCF7L2",              # CD16+mono（非经典）
  "CD68", "CD163", "C1QA",# macro
  "CX3CR1","NR4A3",
  "LILRA4", "IL3RA", "GZMB",                         # pDC
  "FCER1A", "CD1E", "CD1C", "CD1D", "WFDC21P",                  # mDC
  "CLEC9A", "CADM1", "XCR1", "FLT3", "IDO1",                    # mDC1
  "FCER1A","CD1C", "CLEC10A", "CD1E", "HLA-DQA1",               # mDC2
  "LAMP3", "FSCN1", "CCR7",                                     # mDC3
  "EBI3", "LAMP3", "CCR7",                            # mregDC
  "KIT", "TPSAB1", "CPA3",                                      # mast
  "FCER1A", "HDC", "MS4A2", "CLC", "MCPT8",  # basophil
  "FCGR3B", "CSF3R", "CXCR2", "G0S2", "S100A9", "S100A8", "ITGAM", "LTF", "CEACAM8" # neutrophils
)


VlnPlot(object_Myeloid,features = marker_list_myeloid,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2") + NoLegend()
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "group",raster=F)
DimPlot(object_Myeloid, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()

Myeloidjson <- rjson::fromJSON(file = "~/Single Cell/BPH_LS/nameing json/Myeloid cell.json")
object_Myeloid <- label_clusters(object_Myeloid,Myeloidjson)
DimPlot(object = object_Myeloid,reduction ="harmony.secondary.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
VlnPlot(object_Myeloid,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
object_Myeloid<- subset(object_Myeloid, subset = secondary_type != "delete")
saveRDS(object_Myeloid,file = "BPH_LS.secondary.annotation.Myeloid.rds", compress = FALSE)


