library(Seurat)
colnames(object_merged@meta.data)
VlnPlot(object_merged, c("nCount_RNA", "nFeature_RNA","percent.mt"), flip = T, pt.size = 0,group.by = "clusters.harmony.major.pcs40.umap.dims30.res0.8") + NoLegend()
object_merged$barcodes <- rownames(object_merged@meta.data)
cell_quality_data <- as.data.frame(object_merged@meta.data[c("barcodes", "nCount_RNA", "nFeature_RNA", identity)])
# split.by sample
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims40", group.by = "clusters.harmony.major.pcs40.umap.dims40.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims40", group.by = "clusters.harmony.major.pcs40.umap.dims40.res0.5",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims40", group.by = "clusters.harmony.major.pcs40.umap.dims40.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims40", group.by = "clusters.harmony.major.pcs40.umap.dims40.res0.8",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims40", group.by = "clusters.harmony.major.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims40", group.by = "clusters.harmony.major.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()

DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims30", group.by = "clusters.harmony.major.pcs40.umap.dims30.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims30", group.by = "clusters.harmony.major.pcs40.umap.dims30.res0.5",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims30", group.by = "clusters.harmony.major.pcs40.umap.dims30.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims30", group.by = "clusters.harmony.major.pcs40.umap.dims30.res0.8",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims30", group.by = "clusters.harmony.major.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs40.umap.dims30", group.by = "clusters.harmony.major.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()

DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims30", group.by = "clusters.harmony.major.pcs30.umap.dims30.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims30", group.by = "clusters.harmony.major.pcs30.umap.dims30.res0.5",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims30", group.by = "clusters.harmony.major.pcs30.umap.dims30.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims30", group.by = "clusters.harmony.major.pcs30.umap.dims30.res0.8",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims30", group.by = "clusters.harmony.major.pcs30.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims30", group.by = "clusters.harmony.major.pcs30.umap.dims30.res1.2",label = T)+ NoLegend()

DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims20", group.by = "clusters.harmony.major.pcs30.umap.dims20.res0.5", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims20", group.by = "clusters.harmony.major.pcs30.umap.dims20.res0.5",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims20", group.by = "clusters.harmony.major.pcs30.umap.dims20.res0.8", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims20", group.by = "clusters.harmony.major.pcs30.umap.dims20.res0.8",label = T)+ NoLegend()
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims20", group.by = "clusters.harmony.major.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_merged, reduction = "harmony.major.pcs30.umap.dims20", group.by = "clusters.harmony.major.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
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
VlnPlot(object_merged, marker_list, flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.major.pcs40.umap.dims30.res0.8") + NoLegend()
# rename cell
object_merged[["primary_type"]] <- "NA"
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(8,9,15,16), "Epithelia", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(2), "Endothelia", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(7,10,19), "Fibroblast", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(0,3,4,5,6,11), "T/NK cell", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(12), "B cell", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(1,14,17,20), "Myeloid cell", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(18), "Mast cell", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(13), "Unknow cell", object_merged$primary_type )
object_merged$primary_type <- ifelse(object_merged$clusters.harmony.major.pcs40.umap.dims30.res0.8 %in% c(21,22), "Delete", object_merged$primary_type )


table(object_merged@meta.data$primary_type)
object_merged <- subset(object_merged, primary_type!= "Delete" ) 
DimPlot(object_merged,reduction="harmony.major.pcs40.umap.dims30",group.by = "primary_type",label = T,raster = F)+NoLegend()
VlnPlot(object_merged, marker_list, flip = T, pt.size = 0, stack = T, group.by = "primary_type") + NoLegend()

saveRDS(object_merged,file="BPH_LS.major.annotation.pc40umap30res0.8.rds")

Idents(object_merged) <- "primary_type"
object_merged_joined <- JoinLayers(object_merged)
unknow_deg <- FindMarkers(object_merged_joined ,ident.1 = "Unknow cell",logfc.threshold= 0.25, min.pct= 0.1,only.pos= T)
write.csv(unknow_deg ,"unknow_markers.csv")
