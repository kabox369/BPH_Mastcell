table(BPH_SLN.major.annotation.pc40umap40res0.8@meta.data$primary_type)
Idents(object = BPH_SLN.major.annotation.pc40umap40res0.8) <- "primary_type"

object_fibroblast<- subset(BPH_SLN.major.annotation.pc40umap40res0.8,idents=c("fibroblast","myofibroblast"))
object_fibroblast <- NormalizeData(object_fibroblast)
object_fibroblast <- FindVariableFeatures(object_fibroblast)
object_fibroblast <- ScaleData(object_fibroblast)
object_fibroblast <- RunPCA(object_fibroblast)
object_fibroblast <- IntegrateLayers(object = object_fibroblast, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                  verbose = FALSE
)
object_fibroblast <- RunUMAP(object_fibroblast, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_fibroblast <- FindNeighbors(object_fibroblast, reduction = "harmony.pcs30", dims = 1:30)
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.3")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.6")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.9")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")

object_fibroblast <- RunUMAP(object_fibroblast, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_fibroblast <- FindNeighbors(object_fibroblast, reduction = "harmony.pcs30", dims = 1:20)
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.3")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.6")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.9")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_fibroblast <- IntegrateLayers(
  object = object_fibroblast, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_fibroblast <- RunUMAP(object_fibroblast, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_fibroblast <- FindNeighbors(object_fibroblast, reduction = "harmony.pcs40", dims = 1:40)
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.3")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.6")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.9")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_fibroblast <- RunUMAP(object_fibroblast, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_fibroblast <- FindNeighbors(object_fibroblast, reduction = "harmony.pcs40", dims = 1:30)
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.3, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.3")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.6, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.6")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 0.9, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.9")
object_fibroblast <- FindClusters(object_fibroblast, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.3",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res0.9",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims40", group.by = "clusters.harmony.pcs40.umap.dims40.res1.2",label = T)+ NoLegend()


DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()


DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.3",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.6",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs30.umap.dims30.res0.9",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims40.res0.6",label = T)+ NoLegend()


DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.3",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.6",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res0.9",label = T)+ NoLegend()
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2", split.by = "samplename", ncol = 5 , raster=F)
DimPlot(object_fibroblast, reduction = "harmony.pcs30.umap.dims20", group.by = "clusters.harmony.pcs30.umap.dims20.res1.2",label = T)+ NoLegend()
saveRDS(object_fibroblast,file = "bph.subclass.myeloid.rds", compress = FALSE)

# 创建向量，确保每个元素都是明确的
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
  "VWF", "PECAM1",         # endo
  "TPM2", "RGS5",          # SMC
  "DCN", "COL1A1","COL1A2", # Fibroblast
  "COL1A2","TAGLN","ACTA2","POSTN","MMP11","MYL9",#myfibroblast
  "CXCL12", "SOD2","PDGFRA", "IL6", "CXCL14",#inflammatory fibroblast
  "CD74","HLA-DRA","HLA-DRB1","C1QA","C1QB","HLA-DPA1","HLA-DQA", #Antigen-presenting fibroblast
  "PI16","DPP4","LY6C1","COl4A1","HSPG2","CAL15A1"#fibroblast progenitors
)

p5 <- VlnPlot(object_fibroblast,features = marker_list_fibroblast,flip = T, pt.size = 0, stack = T, group.by = "clusters.harmony.pcs40.umap.dims30.res1.2") + NoLegend()
p5.1 <- DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2", split.by = "samplename",ncol = 3,raster=F)
p5.2 <- DimPlot(object_fibroblast, reduction = "harmony.pcs40.umap.dims30", group.by = "clusters.harmony.pcs40.umap.dims30.res1.2",label = T)+ NoLegend()
ggsave(plot = p5,filename = "vinplot.BPH_NSL.fibroblast.pcs40.dims30res1.2.pdf",width = 20,height = 20)
ggsave(plot = p5.1,filename = "samplename.UMAP.BPH_NSL.fibroblast.pcs40.dims30res1.2.pdf",width = 15,height = 30)
ggsave(plot = p5.2,filename = "umap.BPH_NSL.fibroblast.pcs40.dims30res1.2.pdf",width = 7,height = 7)

fibroblastjson <- fromJSON(file = "~/Single Cell/BPH_NSL/secondary annotation/naming json reumap/fibroblast.json")
object_fibroblast <- label_clusters(object_fibroblast,fibroblastjson)
DimPlot(object = object_fibroblast,reduction ="harmony.pcs40.umap.dims30",group.by = "secondary_type" ,label = T )
p <- VlnPlot(object_fibroblast,features = marker_list_fibroblast,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()
ggsave(plot = p,filename = "fibroblast subclass annotationpc40dim30res1.2.pdf",width = 5,height = 20)
object_fibroblast<- subset(object_fibroblast, subset = secondary_type != "delete")
saveRDS(object_fibroblast,file = "BPH_SLN.secondary.annotation.fibroblast.rds", compress = FALSE)
VlnPlot(object_fibroblast,features = marker_list_myeloid,flip = T, pt.size = 0, stack = T, group.by ="secondary_type") + NoLegend()

DimPlot(object_fibroblast,group.by = "secondary_type",split.by = "group",pt.size = 0,label = T)+ NoLegend()


table(object_fibroblast@meta.data$secondary_type)
Idents(object_fibroblast)="secondary_type"#先把Idents改回来
object_fibroblast <- JoinLayers(object_fibroblast)
# 找出每个cluster的标记与所有剩余的细胞相比较，only.pos = T只报告阳性细胞
diff_Fibroblast <- FindAllMarkers(object = object_fibroblast, 
                               only.pos = FALSE, 
                               test.use = "wilcox", 
                               slot = "data", 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25 
)
write.csv(diff_Fibroblast,"diff_Fibroblast.csv")
Fibroblast_degs = FindMarkers( object_fibroblast, 
                           logfc.threshold = 0.25,
                           min.pct = 0.1, # 表达比例的阈值设置，小一点找出更多差异基因
                           only.pos = FALSE,
                           ident.1 = "Myofibroblast_RGS5+",
                           ident.2 = "Fibroblast_SOX+")

write.csv(Fibroblast_degs,"diff_fibroblast_2cluster.csv")

