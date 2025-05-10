Idents(object = object_T_NK) <- "secondary_type"
object_NK_likeT<- subset(object_T_NK,idents=c("CD8+Temra_NK like"))
object_NK_likeT <- NormalizeData(object_NK_likeT)
object_NK_likeT <- FindVariableFeatures(object_NK_likeT)
object_NK_likeT <- ScaleData(object_NK_likeT)
object_NK_likeT <- RunPCA(object_NK_likeT)
object_NK_likeT <- IntegrateLayers(object = object_NK_likeT, method = HarmonyIntegration,
                               orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                               verbose = FALSE
)
object_NK_likeT <- RunUMAP(object_NK_likeT, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_NK_likeT <- FindNeighbors(object_NK_likeT, reduction = "harmony.pcs30", dims = 1:30)
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_NK_likeT <- RunUMAP(object_NK_likeT, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_NK_likeT <- FindNeighbors(object_NK_likeT, reduction = "harmony.pcs30", dims = 1:20)
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_NK_likeT <- IntegrateLayers(
  object = object_NK_likeT, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_NK_likeT <- RunUMAP(object_NK_likeT, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_NK_likeT <- FindNeighbors(object_NK_likeT, reduction = "harmony.pcs40", dims = 1:40)
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_NK_likeT <- RunUMAP(object_NK_likeT, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_NK_likeT <- FindNeighbors(object_NK_likeT, reduction = "harmony.pcs40", dims = 1:30)
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_NK_likeT <- FindClusters(object_NK_likeT, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
marker_list_T_NK <- c(
  "MS4A1", "CD79A",        # B
  "CD14", "FCGR3A",        # macro
  "KIT", "CPA3",           # mast
  "LILRA4", "IL3RA",       # pDC
  "CD1C", "CD1E", "FCER1A", # mDC
  "PPBP",           # mega
  "FCER1A", "HDC",         # basophil
  "LTF","CSF3R",         # neutrophils
  "EPCAM", "ACPP",         # epi
  "VWF", "PECAM1",         # endo
  "ACTA2", "MYH11",        # SMC
  "DCN", "LUM",            # fibro
  "KLK3", "ACPP", "MSMB", "KLK2", "EPCAM", # luminal cell
  "KRT5", "KRT14", "KRT15", "TP63", # basal cell
  "CHGA",#NE
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
  "NR4A1","MYADM","PTGER4", #Trm
  "INFG","NR4A3","CCL4", #Tem
  "RORA","RORC","CCR6","IL22","IL23R","IL17A","IL17F","IL26", #Th17
  "LEF1", "SELL","TCF7","CCR7", #Tn
  "GZMB","GZMH","KLRD1","FCGR3A", #NKT
  "TRDV2","TRGV9","TRGV10", #Tgd
  "MKI67"#Proliferative T
)


p1 <- VlnPlot(object_NK_likeT,features = marker_list_T_NK,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs30.umap.dims30.res1.2") + NoLegend()
p1.1 <- DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename",ncol = 3,raster=F)
p1.2 <- DimPlot(object_NK_likeT, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()

ggsave(plot = p1,filename = "vinplot.NK like T.BPH_NSL.pcs30.dims30res1.2.pdf",width = 20,height = 20)
ggsave(plot = p1.1,filename = "samplename.umap.NK like T.BPH_NSL.pcs30.dims30res1.2.pdf",width = 15,height = 30)
ggsave(plot = p1.2,filename = "umap.NK like T.BPH_NSL.pcs30.dims30res1.2.pdf",width = 7,height = 7)
