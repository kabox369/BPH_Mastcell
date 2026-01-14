unique(object_merged$primary_type)
Idents(object = object_merged) <- "primary_type"
object_Epi<- subset(object_merged,idents = "Epithelia")
object_Epi<- subset(object_Epi,subset = samplename!="BPH_L3")
object_Epi <- NormalizeData(object_Epi)
object_Epi <- FindVariableFeatures(object_Epi)
object_Epi <- ScaleData(object_Epi)
object_Epi <- RunPCA(object_Epi)
object_Epi <- IntegrateLayers(object = object_Epi, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs30", ndim = 30,
                                  verbose = FALSE
)
object_Epi <- RunUMAP(object_Epi, reduction = "harmony.secondary.pcs30", dims = 1:30, reduction.name = "harmony.secondary.pcs30.umap.dims30")
object_Epi <- FindNeighbors(object_Epi, reduction = "harmony.secondary.pcs30", dims = 1:30)
object_Epi <- FindClusters(object_Epi, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3")
object_Epi <- FindClusters(object_Epi, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6")
object_Epi <- FindClusters(object_Epi, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9")
object_Epi <- FindClusters(object_Epi, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2")

object_Epi <- RunUMAP(object_Epi, reduction = "harmony.secondary.pcs30", dims = 1:20, reduction.name = "harmony.secondary.pcs30.umap.dims20")
object_Epi <- FindNeighbors(object_Epi, reduction = "harmony.secondary.pcs30", dims = 1:20)
object_Epi <- FindClusters(object_Epi, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3")
object_Epi <- FindClusters(object_Epi, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6")
object_Epi <- FindClusters(object_Epi, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9")
object_Epi <- FindClusters(object_Epi, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2")

object_Epi <- IntegrateLayers(
  object = object_Epi, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs40", ndim = 40,
  verbose = FALSE
)
object_Epi <- RunUMAP(object_Epi, reduction = "harmony.secondary.pcs40", dims = 1:40, reduction.name = "harmony.secondary.pcs40.umap.dims40")
object_Epi <- FindNeighbors(object_Epi, reduction = "harmony.secondary.pcs40", dims = 1:40)
object_Epi <- FindClusters(object_Epi, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3")
object_Epi <- FindClusters(object_Epi, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6")
object_Epi <- FindClusters(object_Epi, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9")
object_Epi <- FindClusters(object_Epi, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2")

object_Epi <- RunUMAP(object_Epi, reduction = "harmony.secondary.pcs40", dims = 1:30, reduction.name = "harmony.secondary.pcs40.umap.dims30")
object_Epi <- FindNeighbors(object_Epi, reduction = "harmony.secondary.pcs40", dims = 1:30)
object_Epi <- FindClusters(object_Epi, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3")
object_Epi <- FindClusters(object_Epi, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6")
object_Epi <- FindClusters(object_Epi, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9")
object_Epi <- FindClusters(object_Epi, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2")

DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()

marker_list_epi <- c( 
  "CD3D", "CD3E", "CD3G",  #T
  "CD40LG", "CD4",         #CD4
  "CD8A", "CD8B", "GZMK",  # CD8
  "MS4A1", "CD79A",        # B
  "NKG7", "KLRF1",         # NK
  "CD14", "FCGR3A",        # macro
  "KIT", "CPA3",           # mast
  "LILRA4", "IL3RA",       # pDC
  "CD1C", "CD1E", "FCER1A", # mDC
   "PF4",           # mega
  "FCER1A", "HDC",         # basophil
  "LTF","CSF3R",         # neutrophils
  "VWF", "PECAM1",         # endo
  "ACTA2", "MYH11",        # SMC
  "DCN", "LUM",            # fibro
  "KLRF1", "KLRD1",  # NK
  "KLK3", "ACPP", "MSMB", "KLK2", "EPCAM", # luminal cell
  "KRT5", "KRT14", "KRT15", "TP63", # basal cell
  "CHGA","SCG2",#NE
  "SCGB1A1", "SCGB3A1", # club cell
  "KRT13", "KRT5", # hillock cell
  "MKI67","PCNA"
)


VlnPlot(object_Epi,features = marker_list_epi,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6") + NoLegend()
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "group",raster=F)
DimPlot(object_Epi, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
VlnPlot(object_Epi,features = marker_list_epi,flip = T, pt.size = 0, stack = T, group.by = "secondary_type") + NoLegend()

Epijson <- rjson::fromJSON(file = "~/Single Cell/BPH_LS/nameing json/Epithelia.json")
object_Epi <- label_clusters(object_Epi,Epijson)
DimPlot(object = object_Epi,reduction ="harmony.secondary.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
VlnPlot(object_Epi,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
object_Epi<- subset(object_Epi, subset = secondary_type != "delete")
saveRDS(object_Epi,file = "BPH_LS.secondary.annotation.Epithelia.rds", compress = FALSE)


