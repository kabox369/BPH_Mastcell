library(Seurat)
library(ggplot2)
library(rjson)
Idents(object = BPH_SLN.major.annotation.pc40umap40res0.8) <- "primary_type"
object_Bcell<- subset(BPH_SLN.major.annotation.pc40umap40res0.8,idents=c("B cell"))
object_Bcell <- NormalizeData(object_Bcell)
object_Bcell <- FindVariableFeatures(object_Bcell)
object_Bcell <- ScaleData(object_Bcell)
object_Bcell <- RunPCA(object_Bcell)
object_Bcell <- IntegrateLayers(object = object_Bcell, method = HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                verbose = FALSE
)
object_Bcell <- RunUMAP(object_Bcell, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_Bcell <- FindNeighbors(object_Bcell, reduction = "harmony.pcs30", dims = 1:30)
object_Bcell <- FindClusters(object_Bcell, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_Bcell <- FindClusters(object_Bcell, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_Bcell <- RunUMAP(object_Bcell, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_Bcell <- FindNeighbors(object_Bcell, reduction = "harmony.pcs30", dims = 1:20)
object_Bcell <- FindClusters(object_Bcell, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_Bcell <- FindClusters(object_Bcell, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_Bcell <- IntegrateLayers(
  object = object_Bcell, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_Bcell <- RunUMAP(object_Bcell, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_Bcell <- FindNeighbors(object_Bcell, reduction = "harmony.pcs40", dims = 1:40)
object_Bcell <- FindClusters(object_Bcell, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_Bcell <- FindClusters(object_Bcell, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_Bcell <- RunUMAP(object_Bcell, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_Bcell <- FindNeighbors(object_Bcell, reduction = "harmony.pcs40", dims = 1:30)
object_Bcell <- FindClusters(object_Bcell, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_Bcell <- FindClusters(object_Bcell, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_Bcell <- FindClusters(object_Bcell, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_Bcell, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_Bcell,file = "bph.contain.transplant.bcell.unsubclass.rds", compress = FALSE)


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


p6 <- VlnPlot(object_Bcell,features = marker_list_Bcell,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims30.res0.9") + NoLegend()
p6.1 <- DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename",ncol = 3,raster=F)
p6.2 <- DimPlot(object_Bcell, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
ggsave(plot = p6,filename = "vinplot.BPH_NSL.Bcell.pcs40.dims30res0.9.pdf",width = 20,height = 20)
ggsave(plot = p6.1,filename = "samplename.UMAP.BPH_NSL.Bcell.pcs40.dims30res0.9.pdf",width = 15,height = 30)
ggsave(plot = p6.2,filename = "umap.BPH_NSL.Bcell.pcs40.dims30res0.9.pdf",width = 7,height = 7)

bcelljson <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/nameing json/B cell.json")
object_Bcell <- label_clusters(object_Bcell,bcelljson)
DimPlot(object = object_Bcell,reduction ="harmony.pcs40.umap.dims30",group.by = "secondary_type" ,label = T )
p <- VlnPlot(object_Bcell,features = marker_list_Bcell,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
ggsave(plot = p,filename = "B cell subclass annotationpc40dim30res0.9.pdf",width = 5,height = 20)
object_Bcell<- subset(object_Bcell, subset = secondary_type != "delete")
saveRDS(object_Bcell,file = "BPH_SLN.secondary.annotation.B cell.rds", compress = FALSE)


DimPlot(object_T_NK,group.by = "secondary_type",split.by = "group",pt.size = 0,label = T)+ NoLegend()

