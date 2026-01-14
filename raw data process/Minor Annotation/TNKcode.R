library(jsonlite)
Idents(object = object_merged) <- "primary_type"
object_T_NK <- subset(object_merged,idents = "T/NK cell")
object_T_NK <- NormalizeData(object_T_NK)
object_T_NK <- FindVariableFeatures(object_T_NK)
object_T_NK <- ScaleData(object_T_NK)
object_T_NK <- RunPCA(object_T_NK)
object_T_NK <- IntegrateLayers(object = object_T_NK, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony.secondary.pcs30", ndim = 30,
                               verbose = FALSE
)
object_T_NK <- RunUMAP(object_T_NK, reduction = "harmony.secondary.pcs30", dims = 1:30, reduction.name = "harmony.secondary.pcs30.umap.dims30")
object_T_NK <- FindNeighbors(object_T_NK, reduction = "harmony.secondary.pcs30", dims = 1:30)
object_T_NK <- FindClusters(object_T_NK, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9")
object_T_NK <- FindClusters(object_T_NK, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2")

object_T_NK <- RunUMAP(object_T_NK, reduction = "harmony.secondary.pcs30", dims = 1:20, reduction.name = "harmony.secondary.pcs30.umap.dims20")
object_T_NK <- FindNeighbors(object_T_NK, reduction = "harmony.secondary.pcs30", dims = 1:20)
object_T_NK <- FindClusters(object_T_NK, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9")
object_T_NK <- FindClusters(object_T_NK, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2")

object_T_NK <- IntegrateLayers(
  object = object_T_NK, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.secondary.pcs40", ndim = 40,
  verbose = FALSE
)
object_T_NK <- RunUMAP(object_T_NK, reduction = "harmony.secondary.pcs40", dims = 1:40, reduction.name = "harmony.secondary.pcs40.umap.dims40")
object_T_NK <- FindNeighbors(object_T_NK, reduction = "harmony.secondary.pcs40", dims = 1:40)
object_T_NK <- FindClusters(object_T_NK, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9")
object_T_NK <- FindClusters(object_T_NK, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2")

object_T_NK <- RunUMAP(object_T_NK, reduction = "harmony.secondary.pcs40", dims = 1:30, reduction.name = "harmony.secondary.pcs40.umap.dims30")
object_T_NK <- FindNeighbors(object_T_NK, reduction = "harmony.secondary.pcs40", dims = 1:30)
object_T_NK <- FindClusters(object_T_NK, resolution = 0.3, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.6, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6")
object_T_NK <- FindClusters(object_T_NK, resolution = 0.9, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9")
object_T_NK <- FindClusters(object_T_NK, resolution = 1.2, cluster.name = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2")

DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims40", group.by = "clusters.harmony.secondary.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs40.umap.dims30", group.by = "clusters.harmony.secondary.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims20", group.by = "clusters.harmony.secondary.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_T_NK,file = "BPH_NSL.subclass.T_NK.rds", compress = FALSE)


marker_list_T_NK <- c(
  "MS4A1", "CD79A",        # B
  "CD14", "FCGR3A",        # macro
  "KIT", "CPA3",           # mast
  "LILRA4", "IL3RA",       # pDC
  "CD1C", "CD1E", "FCER1A", # mDC
  "PPBP", "PF4",           # mega
  "FCER1A", "HDC",         # basophil
  "LTF","CSF3R",         # neutrophils
  "EPCAM", "ACPP",         # epi
  "VWF", "PECAM1",         # endo
  "ACTA2", "MYH11",        # SMC
  "DCN", "LUM",            # fibro
  "KLK3", "ACPP", "MSMB", "KLK2", "EPCAM", # luminal cell
  "KRT5", "KRT14", "KRT15", "TP63", # basal cell
  "CHGA","SCG2",#NE
  "SCGB1A1", "SCGB3A1", # club cell
  "KRT13", "KRT5", # hillock cell
  "NKG7", "GNLY", "KLRF1", "KLRD1", "MKI67", "NCAM1", "FCGR3A" ,
  "CD3D","CD3E","CD3G",#T cell
  "CD8A", "CD8B", "GZMK",  #CD8
  "CD40LG","CD4",#CD4
  "LAYN", "LAG3", "TIGIT", "PDCD1", "HAVCR2", "CTLA4", "TCF1",  #Tex 
  "ITGAE","XCL1", "XCL2", "MYADM", "CD6",  #Trm
  "TBX21", "PRF1", "CX3CR1", "FGFBP2", "KLRG1", "ASCL2", "FCGR3A",  #Temra
  "LEF1",  "TCF7", "CCR7","SELL",  #Tn
  "CXCR3","CXCR4","CD44",#Tem
  "IL7R", "CD28", "GPR183", #Tcm
  "SLC4A10","TRAV1-2","ZBTB16","RORC", "KLRB1",#MAIT
  "CXCL13","HAVCR2","INFG",#Th1
  "CXCR4","CD69","ANXA1", #Tcm
  "TBX21","PRF1","CXC3CR1","GNLY","NKG7",#Temra
  "TOX2","TOX","IL21","GNG4","BCL6","MAGEH1","BTLA","CXCR5","PDCD1", #Tfh
  "FOXP3", #Treg
  "NR4A1","MYADM","PTGER4","CD8+Trm_ZNF683", #Trm
  "INFG","NR4A3","CCL4", #Tem
  "RORA","RORC","CCR6","IL22","IL23R","IL17A","IL17F","IL26", #Th17
  "LEF1", "SELL","TCF7","CCR7", #Tn
  "GZMB","GZMH","KLRD1","FCGR3A", #NKT
  "TRDV2","TRGV9","TRGV10", #Tgd
  "MKI67"#Proliferative T
)


VlnPlot(object_T_NK,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2") + NoLegend()
VlnPlot(object_T_NK,features = secondaryfeature,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2") + NoLegend()
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2", split.by = "samplename",ncol = 3,raster=F)
DimPlot(object_T_NK, reduction = "harmony.secondary.pcs30.umap.dims30", group.by = "clusters.harmony.secondary.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()

T_NKjson <- rjson::fromJSON(file = "~/Single Cell/BPH_LS/nameing json/T_NK.json")
object_T_NK <- label_clusters(object_T_NK,T_NKjson)
DimPlot(object = object_T_NK,reduction ="harmony.secondary.pcs30.umap.dims30",group.by = "secondary_type" ,label = T )
VlnPlot(object_T_NK,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
ggsave(plot = p,filename = "vinplot by cluster BPH_SLN.secondary annotation pc30dim30res1.2.pdf",width = 15,height = 20)
object_T_NK<- subset(object_T_NK, subset = secondary_type != "delete")
saveRDS(object_T_NK,file = "BPH_LS.secondary.annotation.T_NK.rds", compress = FALSE)



