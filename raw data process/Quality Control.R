library(reticulate)
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
maximum_umi <- 20000
minimum_umi <- 600
minimum_gene <- 600
maximum_gene <- 8000
output_path <- "~/Desktop/out_path/BPH_LS"
object_list <- list()
# 1. define the data path list
path_list <- list(
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_L1/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_L2/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_L3/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_L4/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_L5/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_L6/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_S1/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_S2/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_S3/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_S4/filtered_feature_bc_matrix",
  "~/Desktop/BPH.sc.seq /1.CellRanger/BPH_S5/filtered_feature_bc_matrix"
)
for (path in path_list) {
  print(path)
  samplename <- str_split_i(path,pattern = "/",i = 5) # 以斜杠切分，取对应位置
  object.data <- Read10X(path)
  object <- CreateSeuratObject(counts = object.data, project = "object", min.cells = 3, min.features = 200)
  object@meta.data$samplename  <- samplename
  object[["percent.mt"]] <- PercentageFeatureSet( object, pattern = "^[Mm][Tt]")
  object[["percent.rpls"]]<- PercentageFeatureSet( object, pattern = "^[Rr][Pp][Ll][Ss]")
  hb.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  matched_genes <- match(hb.genes, rownames(object))
  hb_genes <- hb.genes[!is.na(matched_genes)]
  hb_genes <- na.omit(hb_genes)
  object[["percent.hb"]] <- PercentageFeatureSet(object, features = hb_genes)
  # 细胞周期
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  object <- NormalizeData(object)
  object <- CellCycleScoring(object = object,g2m.features = g2m.genes,s.features =s.genes,set.ident = TRUE)
  # 双细胞
  use_condaenv("/opt/anaconda3/envs/alan_ai/bin/python", required = TRUE)
  scr_py <- import("scrublet")
  counts_matrix <- object@assays$RNA$counts
  counts_matrix <- t(counts_matrix)
  scrub <- scr_py$Scrublet(counts_matrix)
  doublets_result <- scrub$scrub_doublets(verbose = FALSE)
  doublets <- as.data.frame(doublets_result)
  colnames(doublets) <- c("doublet_score", "is_doublet")
  rownames(doublets) <- rownames(object@meta.data)
  write.csv(doublets, file.path(output_path, paste(samplename, "_scrublet_doublets_score.csv", sep = "")), row.names = F)
  object <- AddMetaData(object, doublets)
  saveRDS(object, file = file.path(output_path, paste(samplename, "_unfiltered.rds", sep = "")))
  # 质控
  object@meta.data$quality  <- "Pass"
  # doublets
  object[["quality"]] <- ifelse(object@meta.data$is_doublet == "TRUE", "Doublet", "Pass")
  
  # umi
  object[["quality"]] <- ifelse(object@meta.data$nCount_RNA < minimum_umi & object@meta.data$quality == "Pass", "Low_nCount", object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$nCount_RNA < minimum_umi & object@meta.data$quality != "Pass", paste("Low_nCount", object@meta.data$quality, sep = ","), object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$nCount_RNA > maximum_umi & object@meta.data$quality == "Pass", "High_nCount", object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$nCount_RNA > maximum_umi & object@meta.data$quality != "Pass", paste("High_nCount", object@meta.data$quality, sep = ","), object@meta.data$quality)
  
  # gene
  object[["quality"]] <- ifelse(object@meta.data$nFeature_RNA < minimum_gene & object@meta.data$quality == "Pass", "Low_nFeature", object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$nFeature_RNA < minimum_gene & object@meta.data$quality != "Pass" & object@meta.data$quality != "Low_nFeature", paste("Low_nFeature", object@meta.data$quality, sep = ","), object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$nFeature_RNA > maximum_gene & object@meta.data$quality == "Pass", "High_nFeature", object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$nFeature_RNA > maximum_gene & object@meta.data$quality != "Pass" & object@meta.data$quality != "High_nFeature", paste("High_nFeature", object@meta.data$quality, sep = ","), object@meta.data$quality)
  
  # MT
  object[["quality"]] <- ifelse(object@meta.data$percent.mt > 25 & object@meta.data$quality == "Pass", "High_MT", object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$percent.mt > 25 & object@meta.data$quality != "Pass" & object@meta.data$quality != "High_MT", paste("High_MT", object@meta.data$quality, sep = ","), object@meta.data$quality)
  
  # HB
  object[["quality"]] <- ifelse(object@meta.data$percent.hb > 10 & object@meta.data$quality == "Pass", "High_Erythrocyte", object@meta.data$quality)
  object[["quality"]] <- ifelse(object@meta.data$percent.hb > 10 & object@meta.data$quality != "Pass" & object@meta.data$quality != "High_Erythrocyte", paste("High_Erythrocyte", object@meta.data$quality, sep = ","), object@meta.data$quality)
  
  object <- subset(object, subset=quality=="Pass")
  saveRDS(object, file = file.path(output_path, paste(samplename, "_filtered.rds", sep = "")))
  object <-RenameCells(object,add.cell.id = object$samplename )
  object_list[[samplename]] <- object
}
object_merge <- merge(object_list[[1]], object_list[-1], merge.data = TRUE)
saveRDS(object_merge, file = file.path(output_path, "object.merged.rds"),compress = FALSE)
object_merge <- JoinLayers(object_merge)
object_merge[["RNA"]] <- split(object_merge[["RNA"]], f = object_merge$samplename)
#splitLayers by samplename
object_merge <- NormalizeData(object_merge)
object_merge <- FindVariableFeatures(object_merge)
object_merge <- ScaleData(object_merge)
object_merge <- RunPCA(object_merge)
object_merge <- IntegrateLayers(object = object_merge, method = HarmonyIntegration,
                                orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                verbose = FALSE
)
object_merge <- RunUMAP(object_merge, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
object_merge <- FindNeighbors(object_merge, reduction = "harmony.pcs30", dims = 1:30)
object_merge <- FindClusters(object_merge, resolution = 0.5, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.5")
object_merge <- FindClusters(object_merge, resolution = 0.8, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.8")
object_merge <- FindClusters(object_merge, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")

object_merge <- RunUMAP(object_merge, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
object_merge <- FindNeighbors(object_merge, reduction = "harmony.pcs30", dims = 1:20)
object_merge <- FindClusters(object_merge, resolution = 0.5, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.5")
object_merge <- FindClusters(object_merge, resolution = 0.8, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.8")
object_merge <- FindClusters(object_merge, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")

object_merge <- IntegrateLayers(
  object = object_merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
  verbose = FALSE
)
object_merge <- RunUMAP(object_merge, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
object_merge <- FindNeighbors(object_merge, reduction = "harmony.pcs40", dims = 1:40)
object_merge <- FindClusters(object_merge, resolution = 0.5, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.5")
object_merge <- FindClusters(object_merge, resolution = 0.8, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.8")
object_merge <- FindClusters(object_merge, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")

object_merge <- RunUMAP(object_merge, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
object_merge <- FindNeighbors(object_merge, reduction = "harmony.pcs40", dims = 1:30)
object_merge <- FindClusters(object_merge, resolution = 0.5, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.5")
object_merge <- FindClusters(object_merge, resolution = 0.8, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.8")
object_merge <- FindClusters(object_merge, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")

object_merge$group <- "NA"
object_merge$group <- ifelse(object_merge$samplename %in% c("BPH_L1", "BPH_L2", "BPH_L3", "BPH_L4", "BPH_L5", "BPH_L6"),
                             "L", "S")
DimPlot(object = object_merge,split.by = "group")
saveRDS(object_merge,file = "bph.rawdata.flitered.rds", compress = FALSE)


