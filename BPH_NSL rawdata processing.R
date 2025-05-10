  library(Seurat)
  object_SLN <- merge(x = LS_object.merged,y = NL_object.merged,project = "object_SLN")
  object_SLN <- JoinLayers(object_SLN)
  object_SLN[["RNA"]] <- split(object_SLN[["RNA"]], f = object_SLN$samplename)
  #splitLayers by samplename
  object_SLN <- NormalizeData(object_SLN)
  object_SLN <- FindVariableFeatures(object_SLN)
  object_SLN <- ScaleData(object_SLN)
  object_SLN <- RunPCA(object_SLN)
  object_SLN <- IntegrateLayers(object = object_SLN, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony.pcs30", ndim = 30,
                                  verbose = FALSE
  )
  object_SLN <- RunUMAP(object_SLN, reduction = "harmony.pcs30", dims = 1:30, reduction.name = "harmony.pcs30.umap.dims30")
  object_SLN <- FindNeighbors(object_SLN, reduction = "harmony.pcs30", dims = 1:30)
  object_SLN <- FindClusters(object_SLN, resolution = 0.5, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.5")
  object_SLN <- FindClusters(object_SLN, resolution = 0.8, cluster.name = "clusters.harmony.pcs30.umap.dims30.res0.8")
  object_SLN <- FindClusters(object_SLN, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims30.res1.2")
  
  object_SLN <- RunUMAP(object_SLN, reduction = "harmony.pcs30", dims = 1:20, reduction.name = "harmony.pcs30.umap.dims20")
  object_SLN <- FindNeighbors(object_SLN, reduction = "harmony.pcs30", dims = 1:20)
  object_SLN <- FindClusters(object_SLN, resolution = 0.5, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.5")
  object_SLN <- FindClusters(object_SLN, resolution = 0.8, cluster.name = "clusters.harmony.pcs30.umap.dims20.res0.8")
  object_SLN <- FindClusters(object_SLN, resolution = 1.2, cluster.name = "clusters.harmony.pcs30.umap.dims20.res1.2")
  
  object_SLN <- IntegrateLayers(
    object = object_SLN, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony.pcs40", ndim = 40,
    verbose = FALSE
  )
  object_SLN <- RunUMAP(object_SLN, reduction = "harmony.pcs40", dims = 1:40, reduction.name = "harmony.pcs40.umap.dims40")
  object_SLN <- FindNeighbors(object_SLN, reduction = "harmony.pcs40", dims = 1:40)
  object_SLN <- FindClusters(object_SLN, resolution = 0.5, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.5")
  object_SLN <- FindClusters(object_SLN, resolution = 0.8, cluster.name = "clusters.harmony.pcs40.umap.dims40.res0.8")
  object_SLN <- FindClusters(object_SLN, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims40.res1.2")
  
  object_SLN <- RunUMAP(object_SLN, reduction = "harmony.pcs40", dims = 1:30, reduction.name = "harmony.pcs40.umap.dims30")
  object_SLN <- FindNeighbors(object_SLN, reduction = "harmony.pcs40", dims = 1:30)
  object_SLN <- FindClusters(object_SLN, resolution = 0.5, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.5")
  object_SLN <- FindClusters(object_SLN, resolution = 0.8, cluster.name = "clusters.harmony.pcs40.umap.dims30.res0.8")
  object_SLN <- FindClusters(object_SLN, resolution = 1.2, cluster.name = "clusters.harmony.pcs40.umap.dims30.res1.2")
  
  object_SLN$group <- "NA"
  table(object_SLN@meta.data$samplename)
  object_SLN$group <- ifelse(object_SLN$samplename %in% c("BPH_L1", "BPH_L2", "BPH_L3", "BPH_L4", "BPH_L5", "BPH_L6"),
                               "Large", object_SLN$group)
  object_SLN$group <- ifelse(object_SLN$samplename %in% c("BPH_S1", "BPH_S2", "BPH_S3", "BPH_S4", "BPH_S5"),
                             "Small", object_SLN$group)
  object_SLN$group <- ifelse(object_SLN$samplename %in% c("N1", "N2", "N3"),
                             "Normal", object_SLN$group)
  DimPlot(object = object_SLN,split.by = "group")
  saveRDS(object_SLN,file = "BPH_NLS.rawdata.flitered.rds", compress = FALSE)
  