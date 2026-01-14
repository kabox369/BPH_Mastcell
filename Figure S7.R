#Figure S7 A
object.list <- list(Large = BPH_Large_EndoMT_Cellchat, 
                    Small = BPH_Small_EndoMT_Cellchat )
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Fibroblast', targets.use = NULL, stacked =T, do.stat = TRUE,x.rotation = 0,
        comparison = c(1, 2),color.use = c("#2CA02C","#FF7F0E"),)

#Figure S7 B-C
pathways.show3 <- c("ICAM") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(BPH_Large_noEndoMT_Cellchat, signaling = pathways.show3, layout = "circle",
                    top = 0.2,vertex.weight.max = 1)

pathways.show4 <- c("VCAM") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(BPH_Large_noEndoMT_Cellchat, signaling = pathways.show4, layout = "circle",
                    top = 0.2,vertex.weight.max = 1)

#Figure S7 D-F
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)

{
#/Users/alancao/Desktop/GSE183676_RAW/1562
  L1 <- Read10X(data.dir = "/Users/alancao/Desktop/GSE183676_RAW/1562")
L1 <- CreateSeuratObject(counts = L1, project = "L1")
L1

L2 <- Read10X(data.dir = "/Users/alancao/Desktop/GSE183676_RAW/1579")
L2 = L2[["Gene Expression"]]
L2 <- CreateSeuratObject(counts = L2, project = "L2")
L2

L3 <- Read10X(data.dir = "/Users/alancao/Desktop/GSE183676_RAW/1595")
L3 = L3[["Gene Expression"]]
L3 <- CreateSeuratObject(counts = L3, project = "L3")
L3

L4 <- Read10X(data.dir = "/Users/alancao/Desktop/GSE183676_RAW/1628")
L4 = L4[["Gene Expression"]]
L4 <- CreateSeuratObject(counts = L4, project = "L4")
L4

L5 <- Read10X(data.dir = "/Users/alancao/Desktop/GSE183676_RAW/1652")
L5 = L5[["Gene Expression"]]
L5 <- CreateSeuratObject(counts = L5, project = "L5")
L5

L6 <- Read10X_h5("/Users/alancao/Desktop/GSE145928_RAW/GSM4337069_BPH327PrGF_Via_filtered_feature_bc_matrix.h5") 
L6 <- CreateSeuratObject(counts = L6, project = "L6")

L7 <- Read10X_h5("/Users/alancao/Desktop/GSE145928_RAW/GSM4337070_BPH340PrGF_Via_filtered_feature_bc_matrix.h5") 
L7 <- CreateSeuratObject(counts = L7, project = "L7")

L8 <- Read10X_h5("/Users/alancao/Desktop/GSE145928_RAW/GSM4337071_BPH342PrF_Via_filtered_feature_bc_matrix.h5") 
L8 <- CreateSeuratObject(counts = L8, project = "L8")

N1 <- Read10X_h5("/Users/alancao/Desktop/GSE145928_RAW/GSM4337424_D17PrTzF_Via_filtered_feature_bc_matrix.h5") 
N1 <- CreateSeuratObject(counts = N1, project = "N1")

N2 <- Read10X_h5("/Users/alancao/Desktop/GSE145928_RAW/GSM4337426_D27PrTzF_Via_filtered_feature_bc_matrix.h5") 
N2 <- CreateSeuratObject(counts = N2, project = "N2")

N3 <- Read10X_h5("/Users/alancao/Desktop/GSE145928_RAW/GSM4337428_D35PrTzF_Via_filtered_feature_bc_matrix.h5") 
N3 <- CreateSeuratObject(counts = N3, project = "N3")

 BPH_GEO <- merge(L1, y = c(L2,L3,L4,L5,L6,L7,L8,N1,N2,N3), 
             add.cell.ids = c('L1','L2','L3','L4','L5','L6','L7',"L8",
                              'N1','N2','N3'), 
             project = "BPH_GEO")
head(BPH_GEO)
}

 BPH_GEO[["percent.mt"]] <- PercentageFeatureSet( BPH_GEO, pattern = "^MT-",assay = 'RNA')
 VlnPlot(BPH_GEO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
 
 BPH_GEO <- subset( BPH_GEO, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25 & nCount_RNA > 500)
 BPH_GEO@meta.data[["group"]] = substr( BPH_GEO@meta.data[["orig.ident"]],1,1)

 BPH_GEO <- NormalizeData( BPH_GEO)
 BPH_GEO <- FindVariableFeatures( BPH_GEO)
 BPH_GEO <- ScaleData( BPH_GEO)
 BPH_GEO <- RunPCA( BPH_GEO)

 BPH_GEO <- FindNeighbors( BPH_GEO, dims = 1:30, reduction = "pca")
 BPH_GEO <- FindClusters( BPH_GEO, resolution = 1, cluster.name = "unintegrated_clusters")

 BPH_GEO <- RunUMAP( BPH_GEO, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot( BPH_GEO, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

gc()
options(future.globals.maxSize = 1000 * 1024^20)
 BPH_GEO <- IntegrateLayers(
  object =  BPH_GEO, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE)

 BPH_GEO <- IntegrateLayers(
  object =  BPH_GEO, 
  method = RPCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  verbose = FALSE)

future::plan("multisession",workers = 24)
options(future.globals.maxSize = 1000*1024^20)

 BPH_GEO <- IntegrateLayers(
  object =  BPH_GEO, 
  method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE)

saveRDS( BPH_GEO, file = ' BPH_GEO.rds')

#####CCA######
 BPH_GEO <- FindNeighbors( BPH_GEO, reduction = "integrated.cca", dims = 1:30)
 BPH_GEO <- FindClusters( BPH_GEO, resolution = 0.2, cluster.name = "cca_clusters")

 BPH_GEO <- RunUMAP( BPH_GEO, reduction = "integrated.cca", 
               dims = 1:30, 
               reduction.name = "umap.cca")
p1 <- DimPlot(
   BPH_GEO,
  reduction = "umap.cca",
  group.by = c("orig.ident", "cca_clusters"),
  combine = FALSE, label.size = 2
)
#####RPCA######
 BPH_GEO <- FindNeighbors( BPH_GEO, reduction = "integrated.rpca", dims = 1:30)
 BPH_GEO <- FindClusters( BPH_GEO, resolution = 0.8, cluster.name = "rpca_clusters")

 BPH_GEO <- RunUMAP( BPH_GEO, reduction = "integrated.rpca", 
               dims = 1:30, 
               reduction.name = "umap.rpca")

p2 <- DimPlot(
   BPH_GEO,
  reduction = "umap.rpca",
  group.by = c("orig.ident","rpca_clusters"),
  combine = FALSE, label.size = 2
)

####harmony######
 BPH_GEO <- FindNeighbors( BPH_GEO, reduction = "harmony", dims = 1:30)
 BPH_GEO <- FindClusters( BPH_GEO, resolution = 0.2, cluster.name = "harmony_clusters")

 BPH_GEO <- RunUMAP( BPH_GEO, reduction = "harmony", 
               dims = 1:30, 
               reduction.name = "umap.harmony")
p3 <- DimPlot(
   BPH_GEO,
  reduction = "umap.harmony",
  group.by = c("orig.ident", "harmony_clusters"),
  combine = FALSE, label.size = 2
)

library(patchwork)
wrap_plots(c( p2, p3), ncol = 2, byrow = F)

####Marker######
 BPH_GEO <- JoinLayers( BPH_GEO)
 BPH_GEO.markers <- FindAllMarkers( BPH_GEO, only.pos = TRUE, min.pct = 0.25, 
                                   logfc.threshold = 0.25)

saveRDS( BPH_GEO,file = ' BPH_GEO.rds')

FeaturePlot( BPH_GEO,features = c('C3'),reduction = "umap.harmony",
            pt.size = 0.5)

marker = c('KRT7')
{
p1 = FeaturePlot( BPH_GEO,features = marker,reduction = "umap.harmony",
            pt.size = 0.5,raster=FALSE)
p2 = VlnPlot( BPH_GEO,features = marker,pt.size=0)
p1 | p2
}

new.cluster.ids <- c("Basal Cells", 
                     "Luminal Cells" , "Myeloid Cells" , "Basal Cells", "T Cells", "Myofibroblast",
                     "Fibroblast","Endothelial Cells", "Basal Cells", "Basal Cells", "Mast Cells",
                     "B Cells", "Luminal Cells", "Basal Cells")

#???ƴ?cell label??tsne??umapͼ
names(new.cluster.ids) <- levels(BPH_GEO)
levels( BPH_GEO)
 BPH_GEO <- RenameIdents( BPH_GEO,new.cluster.ids)
levels( BPH_GEO)
DimPlot(
   BPH_GEO,
  reduction = "umap.harmony",
  combine = FALSE, label.size = 2
)

saveRDS( BPH_GEO, file = 'BPH_GEO.rds')


levels(BPH_GEO)
BPH_MFEndo = WhichCells(BPH_GEO, idents = c("Endothelial Cells"))
BPH_MFEndo = subset(x=BPH_GEO,cells = BPH_MFEndo)
DimPlot(BPH_MFEndo, reduction = "umap.harmony",pt.size=0.5,label = TRUE)
table(BPH_MFEndo@meta.data[["orig.ident"]])

BPH_GEO_EndMT = WhichCells(BPH_MFEndo,expression = PECAM1 >0.5 & ACTA2>0.5)
BPH_GEO_EndMT = subset(x=BPH_MFEndo,cells = BPH_GEO_EndMT)
table(BPH_GEO_EndMT@meta.data[["orig.ident"]])
DimPlot(BPH_GEO_EndMT, reduction = "umap.harmony",pt.size=0.5,label = TRUE)

saveRDS(BPH_GEO_EndMT, file = 'BPH_GEO_EndMT.rds')

FeaturePlot(BPH_GEO_EndMT,features = c('PECAM1','ACTA2'), cols = c("gray", "coral2"),
            min.cutoff = 0.5, reduction = "umap.harmony",pt.size = 1)




