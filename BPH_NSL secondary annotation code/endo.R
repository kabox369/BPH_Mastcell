object_endo<- subset(BPH_SLN.major.annotation.pc40umap40res0.8,idents=c("endothelia","EndoMT cell"))
object_endo <- NormalizeData(object_endo)
object_endo <- FindVariableFeatures(object_endo)
object_endo <- ScaleData(object_endo)
object_endo <- RunPCA(object_endo)
object_endo <- IntegrateLayers(object = object_endo, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                               verbose = FALSE
)
object_endo <- RunUMAP(object_endo, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_endo <- FindNeighbors(object_endo, reduction = "harmony.pcs30", dims = 1:30)
object_endo <- FindClusters(object_endo, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_endo <- FindClusters(object_endo, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_endo <- FindClusters(object_endo, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_endo <- FindClusters(object_endo, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_endo <- RunUMAP(object_endo, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_endo <- FindNeighbors(object_endo, reduction = "harmony.pcs30", dims = 1:20)
object_endo <- FindClusters(object_endo, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_endo <- FindClusters(object_endo, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_endo <- FindClusters(object_endo, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_endo <- FindClusters(object_endo, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_endo <- IntegrateLayers(
  object = object_endo, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_endo <- RunUMAP(object_endo, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_endo <- FindNeighbors(object_endo, reduction = "harmony.pcs40", dims = 1:40)
object_endo <- FindClusters(object_endo, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_endo <- FindClusters(object_endo, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_endo <- FindClusters(object_endo, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_endo <- FindClusters(object_endo, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_endo <- RunUMAP(object_endo, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_endo <- FindNeighbors(object_endo, reduction = "harmony.pcs40", dims = 1:30)
object_endo <- FindClusters(object_endo, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_endo <- FindClusters(object_endo, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_endo <- FindClusters(object_endo, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_endo <- FindClusters(object_endo, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_endo, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_endo,file = "bph.contain.transplant.unsubclass.endo.rds", compress = FALSE)

marker_list_endo <- c(
  "CD3D", "CD3E", "CD3G", # T
  "CD40LG", "CD4",         #CD4
  "CD8A", "CD8B", "GZMK",  # CD8
  "MS4A1", "CD79A",        # B
  "NKG7", "KLRF1",         # NK
  "CD14", "FCGR3A",        # macro
  "KIT", "CPA3",           # mast
  "IL3RA",       # pDC
  "CD1C", "CD1E", "FCER1A", # mDC
  "PF4",           # mega
  "FCER1A", "HDC",         # basophil
  "LTF","CSF3R",         # neutrophils
  "EPCAM", "ACPP",         # epi
  "VWF", "SELE","PECAM1","IFI27", #endo
  "ACTA2","MYH11","MYL9","TPM2", "RGS5",#muscle
  "DCN","COL1A1", "COL1A2","LUM", #fibro
  "GJA5","FBLN5","GJA4", #artery
  "ACKR1","CLU","SELP",#vein
  "CA4","CD36", "RGCC",#cap
  "PROX1", "LYVE1","CCL21", #lymph
  "COL4A1","KDR","ESM1"#tip cell
)

p9 <- VlnPlot(object_endo,features = marker_list_endo,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims40.res0.6") + NoLegend()
p9.1 <- DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename",ncol  = 3,raster=F)
p9.2 <- DimPlot(object_endo, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()

ggsave(plot = p9,filename = "vinplot.BPH_NSL.endo.pcs40.dims40res0.6.pdf",width = 5,height = 20)
ggsave(plot = p9.1,filename = "samplename.UMAP.BPH_NSL.endo.pcs40.dims40res0.6.pdf",width = 15,height = 30)
ggsave(plot = p9.2,filename = "umap.BPH_NSL.endo.pcs40.dims40res0.6.pdf",width = 7,height = 7)

endo_json <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/nameing json/endo.json")
object_endo <- label_clusters(object_endo,endo_json)
DimPlot(object = object_endo,reduction ="harmony.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
VlnPlot(object_endo,features = marker_list_endo,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
object_endo<- subset(object_endo, subset = secondary_type != "delete")
saveRDS(object_endo,file = "BPH_SLN.secondary.annotation.endo.rds", compress = FALSE)


