library(Seurat)
library(ggplot2)
Idents(object = BPH_SLN.major.annotation.pc40umap40res0.8) <- "primary_type"
object_mast<- subset(BPH_SLN.major.annotation.pc40umap40res0.8,
                     subset = primary_type %in% c("mast cell") &samplename!= "BPH_S4")
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

DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_mast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_mast,file = "bph.contain.transplant.bcell.unsubclass.rds", compress = FALSE)


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

VlnPlot(object_mast,features = c("TGFB1","CDKN1A","VEGFA","TPSAB1"),group.by = "group")
ggsave(plot = p3,filename = "vinplot.BPH_NSL.mast cell.pcs40.dims40res0.6.pdf",width = 20,height = 20)

p3 <- VlnPlot(object_mast,features = marker_list_mast,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims40.res0.6") + NoLegend()
p3.1 <- DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename",ncol = 3,raster=F)
p3.2 <- DimPlot(object_mast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
ggsave(plot = p3,filename = "vinplot.BPH_NSL.mast cell.pcs40.dims40res0.6.pdf",width = 20,height = 20)
ggsave(plot = p3.1,filename = "samplename.UMAP.BPH_NSL.mast cell.pcs40.dims40res0.6.pdf",width = 15,height = 30)
ggsave(plot = p3.2,filename = "umap.BPH_NSL.mast cell.pcs40.dims40res0.6.pdf",width = 7,height = 7)

mastcelljson <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/nameing json/mast cell.json")
object_mast <- label_clusters(object_mast,mastcelljson)
DimPlot(object = object_mast,reduction ="harmony.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
p <- VlnPlot(object_mast,features = marker_list_Bcell,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
ggsave(plot = p,filename = "B cell subclass annotationpc40dim30res0.9.pdf",width = 5,height = 20)
object_mast<- subset(object_mast, subset = secondary_type != "delete")
saveRDS(object_mast,file = "BPH_SLN.secondary.annotation.mast cell.rds", compress = FALSE)


DimPlot(object_T_NK,group.by = "secondary_type",split.by = "group",pt.size = 0,label = T)+ NoLegend()
