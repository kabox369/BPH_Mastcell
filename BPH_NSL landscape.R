#major annotation figure
library(ggplot2)
library(Seurat)
library(patchwork)
library(extrafont)
source("~/Single Cell/configs/color.R")
source("~/Single Cell/configs/selfFunction.R")
# 转换元数据中的列为数值型
seurat_object@meta.data$nCount_RNA <- as.numeric(seurat_object@meta.data$nCount_RNA)
seurat_object@meta.data$nFeature_RNA <- as.numeric(seurat_object@meta.data$nFeature_RNA)
seurat_object@meta.data$percent.mt <- as.numeric(seurat_object@meta.data$percent.mt)
unique(seurat_object@meta.data$group)

# 绘制 UMAP 图
p1 <- DimPlot(
  seurat_object, 
  reduction = "UMAP", 
  group.by = "primary_type", 
  label = FALSE, 
  cols = my_cols_major, 
  alpha = 1, 
  pt.size = 2.5, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE,
)

# 显示图形
print(p1)

#质控质量图
p0 <- VlnPlot(seurat_object, c("nCount_RNA", "nFeature_RNA","percent.mt"), flip = T, pt.size = 0,group.by = "secondary_type",raster = T,ncol = 1) + NoLegend()
  ggsave("QC pic secondary.pdf",p0 ,width = 15,height =30)

p00 <- VlnPlot(seurat_object, c("nCount_RNA", "nFeature_RNA","percent.mt"), flip = T, pt.size = 0,group.by = "primary_type",raster = T,ncol = 1) + NoLegend()
  ggsave("QC pic primary.pdf",p00 ,width = 10,height =20)

#umap图
  #大类umap
  
  p1 <- DimPlot(
    seurat_object, 
    reduction = "UMAP", 
    group.by = "primary_type", 
    label = FALSE, 
    cols = my_cols_major, 
    alpha = 1, 
    pt.size = 2.5, 
    raster.dpi = c(2048, 2048), 
    raster = TRUE
  )
  
  # 调整图形样式
  p1 <- p1 + 
    ggtitle("Celltype") +  # 添加图片标题
    theme(
      legend.spacing.y = unit(1.5, "cm"),  # 调整图例项之间的垂直间距
      legend.key.size = unit(1.5, "lines"),  # 调整图例项的整体大小
      axis.title.x = element_text(size = 10, face = "bold"),  # 横轴标题样式，字体调小
      axis.title.y = element_text(size = 10, face = "bold"),  # 纵轴标题样式，字体调小
      legend.text = element_text(size = 18),  # 图例文字大小
      strip.text = element_text(size = 20, face = "bold"),  # 分面标题样式
      legend.position = "right",  # 将图例放置在右侧
      legend.margin = margin(0, 20, 0, 0),  # 调整图例与图形的距离
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # 设置标题样式
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 5))  # 调整图例项中点的大小
    )
  
  # 保存图像
  ggsave("major annotation figure.pdf", p1, width = 10, height = 8, dpi = 600)
  


# 大类umap by group 
  p2 <- DimPlot(
    seurat_object, 
    reduction = "UMAP", 
    group.by = "primary_type", 
    label = FALSE, 
    cols = my_cols_major, 
    alpha = 1, 
    pt.size = 2.5, 
    raster.dpi = c(2048, 2048), 
    raster = TRUE,
    split.by  = "group"
  )
  
p2 <- p2 + 
  ggtitle("Celltype Splited by Group") +  # 添加图片标题
  theme(
    legend.spacing.y = unit(1.5, "cm"),  # 调整图例项之间的垂直间距
    legend.key.size = unit(1.5, "lines"),  # 调整图例项的整体大小
    axis.title.x = element_text(size = 10, face = "bold"),  # 横轴标题样式，字体调小
    axis.title.y = element_text(size = 10, face = "bold"),  # 纵轴标题样式，字体调小
    legend.text = element_text(size = 18),  # 图例文字大小
    strip.text = element_text(size = 20, face = "bold"),  # 分面标题样式
    legend.position = "right",  # 将图例放置在右侧
    legend.margin = margin(0, 20, 0, 0),  # 调整图例与图形的距离
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # 设置标题样式
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 5))  # 调整图例项中点的大小
  )
# 保存图像
ggsave("major annotation figure by group.pdf", p2, width = 30, height = 15, dpi = 1200)

###亚类umap
p3 <- DimPlot(
  seurat_object, 
  reduction = "UMAP", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_secondary, 
  alpha = 1, 
  pt.size = 2, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE
)

# 调整图形样式
p3 <- p3 + 
  ggtitle("Celltype Subsetted") +  # 添加标题
  theme(
    legend.spacing.y = unit(1, "cm"),  # 调整图例项之间的垂直间距
    legend.text = element_text(size = 14),  # 调整图例文字大小
    axis.title.x = element_text(size = 12, face = "bold"),  # 横轴标题样式
    axis.title.y = element_text(size = 12, face = "bold"),  # 纵轴标题样式
    strip.text = element_text(size = 14, face = "bold"),  # 分面标题样式
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 标题样式
  )

# 保存图像
ggsave("secondary annotation figure.pdf", p3, width = 20, height = 15, dpi = 800)
#####亚类umap by group
p4 <- DimPlot(
  seurat_object, 
  reduction = "UMAP", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_secondary, 
  alpha = 1, 
  pt.size = 2.5, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE,
  split.by = "group"
)

# 调整图形样式
p4 <- p4 + 
  ggtitle("Secondary Type Splited by Group") +  # 添加标题
  theme(
    legend.spacing.y = unit(1, "cm"),  # 调整图例项之间的垂直间距
    legend.text = element_text(size = 12),  # 调整图例文字大小
    axis.title.x = element_text(size = 12, face = "bold"),  # 横轴标题样式
    axis.title.y = element_text(size = 12, face = "bold"),  # 纵轴标题样式
    strip.text = element_text(size = 14, face = "bold"),  # 分面标题样式
    panel.spacing = unit(1, "lines"),  # 调整分面图间距
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 标题样式
  )

# 保存图像
ggsave("secondary annotation figure by group.pdf", p4, width = 30, height = 15, dpi = 600)#####亚类 umap by sample

#####亚类group by samplename#####
p5 <- DimPlot(
  seurat_object, 
  reduction = "UMAP", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_secondary, 
  alpha = 1, 
  pt.size = 2.5, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE,
  split.by = "samplename", 
  ncol = 5
)

# 提取 ggplot 对象并调整分面图
p5 <- p5 + 
  ggtitle("Secondary Type Splited by Samplename") +  # 添加标题
  theme(
    legend.spacing.y = unit(1, "cm"),  # 调整图例项之间的垂直间距
    legend.text = element_text(size = 12),  # 调整图例文字大小
    axis.title.x = element_text(size = 12, face = "bold"),  # 横轴标题样式
    axis.title.y = element_text(size = 12, face = "bold"),  # 纵轴标题样式
    strip.text = element_text(size = 14, face = "bold"),  # 分面标题样式
    panel.spacing = unit(1, "lines"),  # 调整分面图间距
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 标题样式
  ) +
  facet_wrap(~samplename, ncol = 5, scales = "fixed")  # 强制所有分面图共享相同的横轴范围
ggsave("secondary annotation figure by samplename.pdf", p5, width = 30, height = 15, dpi = 600)
#####大类umap by samplename
p6 <- DimPlot(
  seurat_object, 
  reduction = "UMAP", 
  group.by = "primary_type", 
  label = FALSE, 
  cols = my_cols_major, 
  alpha = 1, 
  pt.size = 2.5, 
  raster.dpi = c(2048, 2048), 
  raster = TRUE,
  split.by = "samplename", 
  ncol = 5
)

p6 <- p6 + 
  ggtitle("Primary Type Splited by Samplename") +  # 添加标题
  theme(
    legend.spacing.y = unit(1, "cm"),  # 调整图例项之间的垂直间距
    legend.text = element_text(size = 12),  # 调整图例文字大小
    axis.title.x = element_text(size = 12, face = "bold"),  # 横轴标题样式
    axis.title.y = element_text(size = 12, face = "bold"),  # 纵轴标题样式
    strip.text = element_text(size = 14, face = "bold"),  # 分面标题样式
    panel.spacing = unit(1, "lines"),  # 调整分面图间距
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 标题样式
  ) +
  facet_wrap(~samplename, ncol = 5, scales = "fixed")  # 强制所有分面图共享相同的横轴范围

ggsave("major annotation figure by samplename.pdf", p6, width = 30, height = 15, dpi = 600)

######
unique(seurat_object@meta.data$group)
# 定义更美观的颜色方案
my_cols_group <- c("Large" = "#1F77B4",  # 蓝色
                   "Small" = "#FF7F0E",  # 橙色
                   "Normal" = "#2CA02C") # 绿色

# 绘制 UMAP 图
p7 <- DimPlot(
  seurat_object, 
  reduction = "UMAP", 
  group.by = "group", 
  label = FALSE, 
  cols = my_cols_group,  # 使用优化后的颜色
  alpha = 0.8,           # 调整点的透明度
  pt.size = 1.5,         # 调整点的大小
  raster.dpi = c(2048, 2048), 
  raster = TRUE
)

# 调整图形样式
p7 <- p7 + 
  ggtitle("UMAP by Group") +  # 添加标题
  theme(
    legend.spacing.y = unit(0.5, "cm"),  # 调整图例项之间的垂直间距
    legend.text = element_text(size = 10, face = "bold"),  # 调整图例文字大小
    axis.title.x = element_text(size = 12, face = "bold"),  # 横轴标题样式
    axis.title.y = element_text(size = 12, face = "bold"),  # 纵轴标题样式
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # 标题样式
    panel.background = element_rect(fill = "white"),  # 设置背景为白色
    legend.position = "right",  # 将图例放置在右侧
    legend.key = element_rect(fill = "white"),  # 设置图例键背景为白色
    legend.title = element_blank()  # 移除图例标题
  ) +
  guides(
    colour = guide_legend(override.aes = list(size = 4, alpha = 1))  # 调整图例项中点的大小和透明度
  )

# 保存图像
ggsave("umap by group non annotation.pdf", p7, width = 6, height = 6, dpi = 600)# FeaturePlot 
plot_list <- lapply(featureplotfeature, function(feature) {
  p <- FeaturePlot(
    seurat_object,
    features = feature,
    reduction = "UMAP",
    label = FALSE,          # 不显示聚类标签
    alpha = 1,              # 透明度
    pt.size = 2.5,          # 点大小
    raster = TRUE,          # 启用栅格化
    raster.dpi = c(2048, 2048)  # 高分辨率
  ) +
    labs(title = feature) +  # 显式设置标题为基因名
    theme(
      axis.title = element_blank(),  # 移除坐标轴标题
      axis.text = element_blank(),   # 移除坐标轴文本
      axis.ticks = element_blank(),  # 移除坐标轴刻度
      plot.title = element_text(size = 10, hjust = 0.5)  # 调整标题样式
    )
  
  # 保留颜色条（参考浓度）
  return(p)
})

# 整合子图：每行 5 列，统一颜色条
p8 <- wrap_plots(plot_list, ncol = 5) + 
  plot_layout(guides = "collect")  # 将所有颜色条收集到公共区域

# 显示图形
print(p8)

ggsave(
  filename = "featureplot.pdf",
  plot = p8,
  width = 20, 
  height = 20,
  dpi = 300
)

featureplotfeature <- c("CD3D","CD8A","CD40LG","CD4","SLC4A10",
                        "NKG7","GNLY","KLRF1","NCAM1",
                        "CD79A","MS4A1","JCHAIN","MZB1",
                        "HLA-DRA","CD14","FCGR3A",
                        "KIT","CPA3",
                        "ACPP","KLK2",
                        "KRT5", "KRT14","KRT15", 
                        "SCGB3A1",
                        "VWF","PECAM1",
                        "DCN","COL1A1","ACTA2","MYH11")
p9 <- FeaturePlot(seurat_object, features = featureplotfeature,ncol = 5,reduction = "UMAP",
                  label = F, alpha = 1, pt.size = 2.5 ,raster.dpi = c(2048,2048) ,raster = T)
ggsave("featureplot1.pdf",p9 ,width = 40,height = 40,dpi = 600)



unique(seurat_object@meta.data$secondary_type)
"secondary_type" <- Idents(seurat_object) 
desired_secondary_order <- c("CD4+Tnaive_LEF1", "CD4+Trm_CLNK", "CD4+Tem_NR4A2", "CD4+Treg_FOXP3","CD4+Tex_CTLA4",
                   "CD8+Teff_GZMK","CD8+Temra_GNLY","CD8+Trm_ZNF683","CD8+MAIT_SLCA410","Proliferative T",
                   "CD56+CD16- NK","CD56+CD16+ NK" ,
                   "B naive","B memory","Plasma",
                   "Macrophage_CX3CR1","Macrophage_NR4A3","Macrophage_FCN1","pDC","mDC",
                   "Mast cell_JAK2","Mast cell_CDKN1A","Mast cell_CXCL8",
                   "Luminal epi_ACPP","Luminal epi_CPNE4",
                   "Basal epi_TP63","Basal epi_KRT13","Basal epi_NR4A1","Basal epi_TNFRSF12A",
                   "Club cell_SCGB3A1",
                   "Myofibroblast_RGS5","Fibroblast_FBLN1",
                   "Endothelia","Stromal cell_PECAM1+ACTA2+"
                   )  # 替换为你的细胞类型顺序
seurat_object@meta.data$secondary_type <- factor(
  seurat_object@meta.data$secondary_type,
  levels = desired_secondary_order  # 指定顺序
)
secondaryfeature <- c("CD3D","CD3E","CD3G",
             "CCR7","SELL","LEF1","IL7R","TCF7",
             "CD4","CD40LG","NR4A2","CLNK","FOXP3",
             "CD8A","CD8B","NKG7","GNLY","GZMB","GZMA","GZMK","GZMH",
             "ITGAE","PDCD1","CTLA4","LAG3","HAVCR2","ZNF683","MKI67",
             "SLC4A10","TRAV1-2",
             "NCAM1","FCGR3A","KLRF1","KLRD1",
             "MS4A1","CD79A","CD79B",
             "IGHD","TCL1A","FCER2","IGHM",
             "AIM2","TNFRSF13B",
             "JCHAIN","MZB1",
             "LYZ","CD14","CD68","CD163","C1QA","FCN1","CX3CR1","NR4A3",
             "KIT","CPA3","TPSAB1","GATA2","CDKN1A","CDKN2A","TGFB1","VEGFA","JAK2","CXCL8",
             "VWF","PECAM1",
             "DCN","COL1A1","COL1A2","ACTA2","RGS5","FBLN1",
             "EPCAM","ACPP","KLK2","CPNE4",
             "KRT13","KRT5","KRT14","KRT15","TP63","NR4A1","TNFRSF12A",
             "SCGB3A1"
             )

P <- VlnPlot(seurat_object,features = secondaryfeature,cols = zzm60colors,raster = F,pt.size = 0,flip = T, stack = T,alpha = 0.1,group.by ="secondary_type")+ NoLegend()
ggsave("Vinplot secondary type.pdf",P,width = 20,height = 30)



#Dotplot
DotPlot_Custom <- function(seuratObj,
                           genes, # list (grouped) or vector (ungrouped)
                           group.by,
                           coord_flip = FALSE,
                           scale = TRUE,
                           dot.scale = 4,
                           gene_expr_cutoff = -Inf,
                           gene_pct_cutoff = 0,
                           cell_expr_cutoff = -Inf,
                           cell_pct_cutoff = 0,
                           panel.spacing_distance = 0.5,
                           return_data = FALSE) {
  # Required libraries
  require(dplyr)
  require(Seurat)
  require(RColorBrewer)
  require(ggplot2)
  require(cowplot)
  require(tidyverse)
  require(viridis)  # 使用viridis色标
  
  # Parameters
  col.min <- -2.5
  col.max <- 2.5
  dot.min <- 0
  
  # Reformat genes
  if (is.list(genes)) {
    if (is.null(names(genes))) {
      names(genes) <- 1:length(genes)
    }
    genes <- stack(genes)
    features <- genes$values
    features_group <- genes$ind
  } else {
    features <- genes
    features_group <- NULL
  }
  
  # Check for missing genes
  subset <- features %in% rownames(seuratObj)
  if (sum(!subset) > 0) {
    cat(paste0(paste(features[!subset], collapse = ", "), " is missing in gene list\n"))
  }
  if (length(features) == 0) {
    stop("No intersecting genes, please check gene name format.\n")
  }
  
  # Prepare plot input
  data.features <- FetchData(seuratObj, cells = colnames(seuratObj), vars = features, slot = "data")
  data.features$id <- seuratObj@meta.data[[group.by]]
  
  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, 1:(ncol(data.features) - 1), drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, function(x) sum(x > 0) / length(x))
    list(avg.exp = avg.exp, pct.exp = pct.exp)
  })
  names(data.plot) <- unique(data.features$id)
  data.plot <- lapply(names(data.plot), function(x) {
    data.use <- as.data.frame(data.plot[[x]])
    data.use$features.plot <- rownames(data.use)
    data.use$features.plot_show <- features
    data.use$id <- x
    data.use
  })
  data.plot <- do.call(rbind, data.plot)
  
  avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
    data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
    if (scale) {
      data.use <- scale(data.use)
      MinMax(data.use, min = col.min, max = col.max)
    } else {
      log1p(data.use)
    }
  })
  avg.exp.scaled <- as.vector(t(avg.exp.scaled))
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(data.plot$features.plot, levels = unique(data.plot$features.plot))
  data.plot$features.plot_show <- factor(data.plot$features.plot_show, levels = unique(features))
  
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if (is.factor(seuratObj@meta.data[[group.by]])) {
    data.plot$id <- factor(data.plot$id, levels = levels(seuratObj@meta.data[[group.by]]))
  }
  
  if (!is.null(features_group)) {
    data.plot <- data.plot %>%
      left_join(data.frame(features = features, features_group = features_group), by = c("features.plot_show" = "features")) %>%
      mutate(features_group = factor(features_group, levels = unique(features_group)))
  }
  
  # Filter genes/cells
  filter_cell <- data.plot %>%
    group_by(id) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp)) %>%
    filter(max_exp >= cell_expr_cutoff & max_pct >= cell_pct_cutoff)
  filter_gene <- data.plot %>%
    group_by(features.plot) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp)) %>%
    filter(max_exp >= gene_expr_cutoff & max_pct >= gene_pct_cutoff)
  
  data.plot <- data.plot %>%
    filter(features.plot %in% filter_gene$features.plot) %>%
    filter(id %in% filter_cell$id)
  
  # Plot
  plot <- ggplot(data.plot, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # 使用蓝色到红色的渐变色
    guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      text = element_text(size = 12),
      plot.margin = unit(c(1, 1, 1, 1), "char"),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed", size = 0.2)
    ) +
    labs(x = "", y = "")
  
  if (coord_flip) {
    plot <- plot + coord_flip()
  }
  
  if (!is.null(features_group)) {
    plot <- plot +
      facet_grid(. ~ features_group, scales = "free_x", space = "free_x", switch = "y") +
      theme(panel.spacing = unit(panel.spacing_distance, "lines"), strip.background = element_blank())
  }
  
  if (return_data) {
    return(list(plot = plot, data = data.plot))
  } else {
    return(plot)
  }
}

# 示例使用
genes <- c("EPCAM","KLK2","MSMB","ACPP",
           "KRT5","KRT14","KRT15","TP63",
           "SCGB1A1","SCGB3A1",
           "VWF","SELE","PECAM1",
           "DCN","COL1A1","LUM",
           "ACTA2", "MYL9","RGS5",
           "PTPRC","CD3D","CD3E","CD3G",
           "MS4A1","CD79A","MZB1",
           "NKG7","GNLY","NCAM1", "KLRF1","KLRD1",
           "LYZ","FCGR3A","CD14",
           "KIT","CPA3","TPSAB1")

group.by <- "primary_type" # 分组依据（如细胞类型）

# 绘制点图
p <- DotPlot_Custom(seuratObj = seurat_object, genes = genes, group.by = group.by)
ggsave("dotplot_primary_type.pdf",p,width = 15,height = 10)
