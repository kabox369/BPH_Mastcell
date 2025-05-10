library(Seurat)
library(ggplot2)
table(BPH_SLN.major.annotation.pc40umap40res0.8@meta.data$primary_type)
Idents(object = BPH_SLN.major.annotation.pc40umap40res0.8) <- "primary_type"
object_myeloid<- subset(BPH_SLN.major.annotation.pc40umap40res0.8,idents=c("myeloid cell"))
object_myeloid <- NormalizeData(object_myeloid)
object_myeloid <- FindVariableFeatures(object_myeloid)
object_myeloid <- ScaleData(object_myeloid)
object_myeloid <- RunPCA(object_myeloid)
object_myeloid <- IntegrateLayers(object = object_myeloid, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                  verbose = FALSE
)
object_myeloid <- RunUMAP(object_myeloid, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_myeloid <- FindNeighbors(object_myeloid, reduction = "harmony.pcs30", dims = 1:30)
object_myeloid <- FindClusters(object_myeloid, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_myeloid <- FindClusters(object_myeloid, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")

object_myeloid <- RunUMAP(object_myeloid, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_myeloid <- FindNeighbors(object_myeloid, reduction = "harmony.pcs30", dims = 1:20)
object_myeloid <- FindClusters(object_myeloid, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_myeloid <- FindClusters(object_myeloid, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_myeloid <- IntegrateLayers(
  object = object_myeloid, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_myeloid <- RunUMAP(object_myeloid, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_myeloid <- FindNeighbors(object_myeloid, reduction = "harmony.pcs40", dims = 1:40)
object_myeloid <- FindClusters(object_myeloid, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_myeloid <- FindClusters(object_myeloid, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_myeloid <- RunUMAP(object_myeloid, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_myeloid <- FindNeighbors(object_myeloid, reduction = "harmony.pcs40", dims = 1:30)
object_myeloid <- FindClusters(object_myeloid, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_myeloid <- FindClusters(object_myeloid, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_myeloid <- FindClusters(object_myeloid, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()


DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_myeloid, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_myeloid,file = "bph.subclass.myeloid.rds", compress = FALSE)

# 创建向量，确保每个元素都是明确的
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
  "CD68", "CD163", "C1QA", "LAMP4",                             # macro
  "LILRA4", "IL3RA", "GZMB", "SERPINF4",                        # pDC
  "FCER1A", "CD1E", "CD1C", "CD1D", "WFDC21P",                  # mDC
  "CLEC9A", "CADM1", "XCR1", "FLT3", "IDO1",                    # mDC1
  "FCER1A","CD1C", "CLEC10A", "CD1E", "HLA-DQA1",               # mDC2
  "LAMP3", "FSCN1", "CCR7",                                     # mDC3
  "EBI3", "LAMP3", "CCR7", "CCL12",                             # mregDC
  "KIT", "TPSAB1", "CPA3",                                      # mast
  "FCER1A", "HDC", "MS4A2", "CLC", "MCPT8", "PRSS34", "CD200R3", # basophil
  "FCGR3B", "CSF3R", "CXCR2", "G0S2", "S100A9", "S100A8", "ITGAM", "LTF", "CEACAM8" # neutrophils
)


p6 <- VlnPlot(object_myeloid,features = marker_list_myeloid,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims40.res1.2") + NoLegend()
p6.1 <- DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename",ncol = 3,raster=F)
p6.2 <- DimPlot(object_myeloid, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()
ggsave(plot = p6,filename = "vinplot.BPH_NSL.myeloid.pcs40.dims40res1.2.pdf",width = 20,height = 20)
ggsave(plot = p6.1,filename = "samplename.UMAP.BPH_NSL.myeloid.pcs40.dims40res1.2.pdf",width = 15,height = 30)
ggsave(plot = p6.2,filename = "umap.BPH_NSL.myeloid.pcs40.dims40res1.2.pdf",width = 7,height = 7)

myeloidjson <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/nameing json/myeloid.json")
object_myeloid <- label_clusters(object_myeloid,myeloidjson)
DimPlot(object = object_myeloid,reduction ="harmony.pcs40.umap.dims40",group.by = "secondary_type" ,label = T )
VlnPlot(object_myeloid,features = marker_list_myeloid,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
ggsave(plot = p,filename = "B cell subclass annotationpc40dim30res0.9.pdf",width = 5,height = 20)
object_myeloid<- subset(object_myeloid, subset = secondary_type != "delete")
saveRDS(object_myeloid,file = "BPH_SLN.secondary.annotation.myeloidcell_final.rds", compress = FALSE)
VlnPlot(object_myeloid,features = marker_list_myeloid,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()

DimPlot(object_myeloid,group.by = "secondary_type",split.by = "group",pt.size = 0,label = T)+ NoLegend()
