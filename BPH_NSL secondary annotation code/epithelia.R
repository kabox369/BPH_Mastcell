Idents(object = BPH_SLN.major.annotation.pc40umap40res0.8) <- "primary_type"
table(BPH_SLN.major.annotation.pc40umap40res0.8@meta.data$primary_type)
object_epi <- subset(BPH_SLN.major.annotation.pc40umap40res0.8, 
                        subset = primary_type %in% c("basal epi", "club cell", "luminal epi") &samplename!= "BPH_L3")
object_epi <- NormalizeData(object_epi)
object_epi <- FindVariableFeatures(object_epi)
object_epi <- ScaleData(object_epi)
object_epi <- RunPCA(object_epi)
object_epi <- IntegrateLayers(object = object_epi, method = HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                verbose = FALSE
)
object_epi <- RunUMAP(object_epi, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_epi <- FindNeighbors(object_epi, reduction = "harmony.pcs30", dims = 1:30)
object_epi <- FindClusters(object_epi, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_epi <- FindClusters(object_epi, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_epi <- FindClusters(object_epi, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_epi <- FindClusters(object_epi, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_epi <- RunUMAP(object_epi, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_epi <- FindNeighbors(object_epi, reduction = "harmony.pcs30", dims = 1:20)
object_epi <- FindClusters(object_epi, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_epi <- FindClusters(object_epi, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_epi <- FindClusters(object_epi, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_epi <- FindClusters(object_epi, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_epi <- IntegrateLayers(
  object = object_epi, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_epi <- RunUMAP(object_epi, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_epi <- FindNeighbors(object_epi, reduction = "harmony.pcs40", dims = 1:40)
object_epi <- FindClusters(object_epi, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_epi <- FindClusters(object_epi, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_epi <- FindClusters(object_epi, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_epi <- FindClusters(object_epi, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_epi <- RunUMAP(object_epi, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_epi <- FindNeighbors(object_epi, reduction = "harmony.pcs40", dims = 1:30)
object_epi <- FindClusters(object_epi, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_epi <- FindClusters(object_epi, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_epi <- FindClusters(object_epi, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_epi <- FindClusters(object_epi, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_epi, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_epi,file = "bph.contain.transplant.bcell.unsubclass.rds", compress = FALSE)


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
  "PPBP", "PF4",           # mega
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
  "KRT13", "KRT5" # hillock cell
)


p4 <- VlnPlot(object_epi,features = marker_list_epi,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims30.res0.9") + NoLegend()
p4.1 <- DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename",ncol = 3,raster=F)
p4.2 <- DimPlot(object_epi, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
ggsave(plot = p4,filename = "vinplot.BPH_NSL.epi.pcs40.dims30res0.9.pdf",width = 20,height = 20)
ggsave(plot = p4.1,filename = "samplename.UMAP.BPH_NSL.epi.pcs40.dims30res0.9.pdf",width = 15,height = 30)
ggsave(plot = p4.2,filename = "umap.BPH_NSL.epi.pcs40.dims30res0.9.pdf",width = 7,height = 7)

epijson <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/nameing json/epi.json")
object_epi <- label_clusters(object_epi,epijson)
DimPlot(object = object_epi,reduction ="harmony.pcs40.umap.dims30",group.by = "secondary_type" ,label = T )
p <- VlnPlot(object_epi,features = marker_list_Bcell,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
ggsave(plot = p,filename = "B cell subclass annotationpc40dim30res0.9.pdf",width = 5,height = 20)
object_epi<- subset(object_epi, subset = secondary_type != "delete")
saveRDS(object_epi,file = "BPH_SLN.secondary.annotation.epi.rds", compress = FALSE)
DimPlot(object_epi,group.by = "secondary_type",split.by = "group",pt.size = 0,label = T)+ NoLegend()

Idents(object_epi)="secondary_type"
object_epi <- JoinLayers(object_epi)
luminal_degs = FindMarkers( object_epi, 
                               logfc.threshold = 0.25,
                               min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                               only.pos = FALSE,
                               ident.1 = "Luminal epi",
                               ident.2 = "Luminal epi_EPCAM-")

write.csv(luminal_degs,"diff_epi_2cluster2.csv")
