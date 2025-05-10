library(Seurat)
colnames(object_SLN@meta.data)
VlnPlot(object_SLN, c("nCount_RNA", "nFeature_RNA","percent.mt"), flip = T, pt.size = 0,group.by = "clusters.harmony.pcs40.umap.dims30.res0.8") + NoLegend()
object_SLN$barcodes <- rownames(object_SLN@meta.data)
cell_quality_data <- as.data.frame(object_SLN@meta.data[c("barcodes", "nCount_RNA", "nFeature_RNA", identity)])
# split.by sample
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.5",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.8",label = T,cols = zzmcolors)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()

DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.5",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.8",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()

DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.5",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.8",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()

DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.5",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.8",label = T)+ NoLegend()
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_SLN, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
# major annotation 
marker_list <- c(
  "PTPRC", #immune cell
  "CD3D", "CD3E", "CD3G",#T cell
  "CD8A","GZMK","CD8B","CCL4L2","ITM2C",#CD8T
  "CD4","CD40LG",#CD4T
  "FOXP3",#Treg
  "MS4A1", "CD79A","CD79B","CD19","MME","MZB1","JCHAIN",#B cell
  "NKG7","GNLY","KLRD1","KLRF1","MKI67", "NCAM1","FCGR3A",#NK
  "MX1","CST3","CSFR1","C1QA","FCGR3A","CD14","LYZ",#macro/mono
  "KIT","CPA3","MS4A1","TPSAB1","TPSB2","CDK15","GATA2",#mast cell
  "CD1E","CD1C","WFDC21P","FCER1A",#pDC
  "FCER1A","CD1E","CD1C","CD1D","WFDC21P",#mDC
  "PPBP", "PF4",#Mega
  "FCGR3B","CSF3R","CXCR2","G0S2","S100A9","S100A8", "ITGAM", "LTF", "CEACAM8",#neutrophils
  "FCER1A","HDC", "MS4A2","CLC", "MCPT8", "PRSS34","CD200R3",#basophil
  "KLK3", "ACPP", "MSMB", "KLK2", "EPCAM", # luminal cell
  "KRT5", "KRT14", "KRT15", "TP63", # basal cell
  "CHGA","SCG2",#NE
  "SCGB1A1", "SCGB3A1", # club cell
  "KRT13", "KRT5", # hillock cell
  "PECAM1", "VWF", "SELE", "IFI27", # Endo
  "ACTA2", "MYH11", "MYL9", "TPM2", "RGS5", # SMC
  "DCN", "COL1A1","COL1A2" # Fibro
)
VlnPlot(object_SLN, marker_list, flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims40.res0.8") + NoLegend()
# rename cell
object_SLN[["primary_type"]] <- "NA"
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(7,16), "luminal epi", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(0,5,11,19), "basal epi", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(21), "club cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(8,9), "endothelia", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(4), "fibroblast", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(6), "myofibroblast", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(20), "EndoMT cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(1,3,12,22,24), "T cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(15), "NK", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(10), "NKT cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(14), "B cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(2,13), "myeloid cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(18), "mast cell", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(17), "unknow", object_SLN$primary_type )
object_SLN$primary_type <- ifelse(object_SLN$clusters.harmony.pcs40.umap.dims40.res0.8 %in% c(23), "delete", object_SLN$primary_type )


table(object_SLN@meta.data$primary_type)
object_SLN <- subset(object_SLN, primary_type!= "delete" ) 
DimPlot(object_SLN,reduction="harmony.pcs40.umap.dims40",group.by = "primary_type",label = T ,cols =  )+NoLegend()
VlnPlot(object_SLN, marker_list, flip = T, pt.size = 0, stack = T, group.by = "primary_type") + NoLegend()

saveRDS(object_SLN,file="BPH_SLN.major.annotation.pc40umap40res0.8.rds")




