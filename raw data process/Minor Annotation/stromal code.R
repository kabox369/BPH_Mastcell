object_Fibro <-
unique(object_merged$primary_type)
Idents(object = object_merged) <- "primary_type"
object_Fibro <- subset(object_merged,idents = c("Fibroblast"))
object_Fibro <- NormalizeData(object_Fibro)
object_Fibro <- FindVariableFeatures(object_Fibro)
object_Fibro <- ScaleData(object_Fibro)
object_Fibro <- RunPCA(object_Fibro)
object_Fibro <- IntegrateLayers(object = object_Fibro, method = HarmonyIntegration,
                              orig.reduction = "pca", new.reduction = "harmony.secondary.pcs30", ndim = 30,
                              verbose = FALSE
)
object_Fibro <- RunUMAP(object_Fibro, reduction = "harmony.secondary.pcs30", dims = 1:30, reduction.name = "harmony.secondary.pcs30.umap.dims30")
object_Fibro <- FindNeighbors(object_Fibro, reduction = "harmony.secondary.pcs30", dims = 1:30)
object_Fibro <- FindClusters(object_Fibro, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9")
object_Fibro <- FindClusters(object_Fibro, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2")

object_Fibro <- RunUMAP(object_Fibro, reduction = "harmony.secondary.pcs30", dims = 1:20, reduction.name = "harmony.secondary.pcs30.umap.dims20")
object_Fibro <- FindNeighbors(object_Fibro, reduction = "harmony.secondary.pcs30", dims = 1:20)
object_Fibro <- FindClusters(object_Fibro, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9")
object_Fibro <- FindClusters(object_Fibro, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2")

object_Fibro <- IntegrateLayers(
  object = object_Fibro, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs40", ndim = 40,
  verbose = FALSE
)
object_Fibro <- RunUMAP(object_Fibro, reduction = "harmony.secondary.pcs40", dims = 1:40, reduction.name = "harmony.secondary.pcs40.umap.dims40")
object_Fibro <- FindNeighbors(object_Fibro, reduction = "harmony.secondary.pcs40", dims = 1:40)
object_Fibro <- FindClusters(object_Fibro, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9")
object_Fibro <- FindClusters(object_Fibro, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2")

object_Fibro <- RunUMAP(object_Fibro, reduction = "harmony.secondary.pcs40", dims = 1:30, reduction.name = "harmony.secondary.pcs40.umap.dims30")
object_Fibro <- FindNeighbors(object_Fibro, reduction = "harmony.secondary.pcs40", dims = 1:30)
object_Fibro <- FindClusters(object_Fibro, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6")
object_Fibro <- FindClusters(object_Fibro, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9")
object_Fibro <- FindClusters(object_Fibro, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2")

DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()

marker_list_fibroblast <- c(
  "CD3D", "CD3E", "CD3G",  #T
  "CD40LG", "CD4",         #CD4
  "CD8A", "CD8B", "GZMK",  # CD8
  "MS4A1", "CD79A",        # B
  "NKG7", "KLRF1","KLRD1", # NK
  "CD14", "FCGR3A",        # mono
  "CD68","CD163",          #macro
  "KIT", "CPA3","TPSAB1",  # mast
  "LILRA4", "IL3RA",       # pDC
  "CD1C", "CD1E", "FCER1A", # mDC
  "PPBP",                  # mega
  "FCER1A", "HDC",         # basophil
  "LTF","CSF3R",           # neutrophils
  "EPCAM", "ACPP",         # epi
  "VWF", "SELE","PECAM1","IFI27", #endo
  "ACTA2","MYH11","MYL9","TPM2", "RGS5",#muscle
  "DCN", "COL1A1","COL1A2", # Fibroblast
  "COL1A2","TAGLN","ACTA2","POSTN","MMP11","MYL9",#myfibroblast
  "CXCL12", "SOD2","PDGFRA", "IL6", "CXCL14",#inflammatory fibroblast
  "CD74","HLA-DRA","HLA-DRB1","C1QA","C1QB","HLA-DPA1","HLA-DQA", #Antigen-presenting fibroblast
  "PI16","DPP4","LY6C1","COl4A1","HSPG2","CAL15A1"#fibroblast progenitors
)

VlnPlot(object_Fibro,features = marker_list_fibroblast,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6") + NoLegend()
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "group",raster=F)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "secondary_type", split.by = "group",raster=F,ncol = 4)
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
VlnPlot(object_Fibro,features = marker_list_epi,flip = T, pt.size = 0, stack = T, group.by = "secondary_type") + NoLegend()

fibrojson <- rjson::fromJSON(file = "~/Single Cell/BPH_LS/nameing json/fibroblast.json")
object_Fibro <- label_clusters(object_Fibro,fibrojson)
DimPlot(object = object_Fibro,reduction ="harmony.secondary.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
VlnPlot(object_Fibro,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
object_Fibro<- subset(object_Fibro, subset = secondary_type != "delete")
DimPlot(object_Fibro, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "secondary_type", label = F, 
        cols = my_cols_secondary2, raster = F)
saveRDS(object_Fibro,file = "BPH_LS.secondary.annotation.Fibroblast.rds", compress = FALSE)


