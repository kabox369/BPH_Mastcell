library(dplyr)
library(assertthat)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
source("~/Single Cell/configs/color.R")
source("~/Single Cell/configs/selfFunction.R")

# 假设 `meta` 是你的数据框，包含了所有的 meta 信息
meta <- seurat_object@meta.data
sample_list <- unique(meta$samplename)
celltype_list <- unique(meta$primary_type)
group_list <- unique(meta$group)

# 定义函数
cell_count <- function(meta_data) {
  # 初始化结果框
  freq_result <- data.frame(matrix(0, nrow = length(celltype_list), ncol = 0))
  rownames(freq_result) <- celltype_list
  
  # 遍历每个分组
  for (group in group_list) {
    message("Processing group: ", group)
    
    # 筛选当前分组的样本
    sample_group <- meta_data[meta_data$group == group, ]
    samples <- unique(sample_group$samplename)
    
    # 遍历每个样本
    for (sample in samples) {
      message("Processing sample: ", sample)
      
      # 筛选当前样本的数据
      person_data <- sample_group[sample_group$samplename == sample, ]
      
      # 计算每种细胞类型的数目
      proportion <- as.data.frame(table(person_data$primary_type)) %>%
        column_to_rownames("Var1")
      group_name <- str_c(group, sample, sep = '_')
      colnames(proportion) <- group_name
      
      # 如果细胞类型数量不足，补齐缺失的细胞类型
      if (nrow(proportion) != length(celltype_list)) {
        message('Adding missing cell types with 0 counts...')
        missing_celltypes <- setdiff(celltype_list, rownames(proportion))
        add_frame <- data.frame(matrix(0, nrow = length(missing_celltypes), ncol = 1))
        rownames(add_frame) <- missing_celltypes
        colnames(add_frame) <- group_name
        proportion <- rbind(proportion, add_frame)
      }
      
      # 确保细胞类型顺序一致
      proportion <- proportion[celltype_list, , drop = FALSE]
      colnames(proportion) <- group_name
      
      # 合并到结果框
      freq_result <- cbind(freq_result, proportion)
    }
  }
  
  # 结果校验：统计meta里面的细胞总数量以及分类计数后的细胞总数量是否匹配
  freq_cache <- rowSums(freq_result)
  raw_freq <- table(meta_data$primary_type)
  
  for (celltype in names(freq_cache)) {
    assert_that(
      freq_cache[[celltype]] == raw_freq[[celltype]],
      msg = paste(celltype, "cell numbers do not match!")
    )
  }
  
  return(freq_result)
}

# 调用函数
cellcount <- cell_count(meta)

# 保存结果
write.csv(cellcount, file = "cell_counts_primary.csv")

## 计算比例
proportion <- apply(cellcount, 2, function(x) x / sum(x))
proportion <- as.data.frame(t(proportion))


proportion <- proportion %>% 
  dplyr::mutate(Group = str_split_fixed(rownames(proportion), '_', n=2)[,1]) %>% 
  dplyr::mutate(Name = str_split_fixed(rownames(proportion), '_', n=2)[,2])


proportion_plot <- gather(proportion, Celltype, Value, -Name, -Group)
proportion_plot$Value <- round(proportion_plot$Value * 100, 2)

proportion_plot$Group <- factor(proportion_plot$Group, levels = c("Large", "Small", "Normal"))

write.csv(proportion, file = "cell_proportion_primary.csv")

unique(proportion_plot$Celltype)

desired_secondary_order <- c(
  "CD4+Tnaive_LEF1", "CD4+Trm_CLNK", "CD4+Tem_NR4A2", "CD4+Treg_FOXP3", "CD4+Tex_CTLA4",
  "CD8+Teff_GZMK", "CD8+Temra_GNLY", "CD8+Trm_ZNF683", "CD8+MAIT_SLCA410", "Proliferative T",
  "CD56+CD16- NK", "CD56+CD16+ NK",
  "B naive", "B memory", "Plasma",
  "Macrophage_CX3CR1", "Macrophage_NR4A3", "Macrophage_FCN1", "pDC", "mDC",
  "Mast cell_JAK2", "Mast cell_CDKN1A", "Mast cell_CXCL8",  # 注意大小写
  "Endothelia", "Stromal cell_PECAM1+ACTA2+",
  "Myofibroblast_RGS5", "Fibroblast_FBLN1",
  "Luminal epi_ACPP", "Luminal epi_CPNE4",
  "Basal epi_TP63", "Basal epi_KRT13", "Basal epi_NR4A1", "Basal epi_TNFRSF12A",
  "Club cell_SCGB3A1"
)
group_order <- c("Large", "Small", "Normal")  # 自定义Group顺序

desired_primary_order <- c("T cell","NK","B cell","Myeloid cell","Mast cell","Epithelia","Endothelia","Fibroblast")




################## Stack figure by samplename###################
legend_order <- desired_primary_order  # 使用 celltype_list 定义的顺序（或自定义）

# 转换因子顺序
proportion_plot$Group <- factor(
  proportion_plot$Group,
  levels = group_order  # 强制分组顺序
)

proportion_plot$Celltype <- factor(
  proportion_plot$Celltype,
  levels = legend_order  # 强制图例顺序
)

# 确保颜色顺序与 legend_order 一致
my_cols_primary <- my_cols_major[legend_order]

# 绘制堆积图
p <- ggplot(proportion_plot, aes(x = Name, y = Value, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +  # 调整柱子宽度
  facet_grid(
    . ~ Group, 
    scales = "free_x", 
    space = "free_x", 
    labeller = labeller(Group = label_value)
  ) +
  scale_fill_manual(
    values = my_cols_primary,
    breaks = legend_order,  # 显式指定图例顺序
    drop = FALSE            # 显示所有因子水平（即使数据中不存在）
  ) +
  labs(
    title = "Cell Type Composition by Sample",
    x = "Samplename",
    y = "Percentage (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.8, "cm"),  # 增大图例项大小
    legend.text = element_text(size = 12),  # 增大图例文字大小
    legend.title = element_text(size = 14),  # 增大图例标题大小
    panel.spacing = unit(1, "lines")  # 调整分面图间距
  )

# 保存图像
ggsave(
  filename = "cell_type_composition_primary_type_Stacked_Bar_Chart.pdf",
  plot = p,
  width = 15,
  height = 12,
  dpi = 300
)
##################### 箱线图模版 proportion figure  by    samplename   ########
p <- ggplot()+ 
  ylab(NULL) +
  theme(panel.background = element_rect(I(0))) +
  theme(axis.line = element_line(color = "black",size = 0.15)) +
  theme(axis.text.x = element_text(size = 7,  color = "black", vjust = 0.5, hjust = 0.5, angle = 90),
        axis.text.y = element_text(size = 7,  color = "black",  vjust = 0.5, hjust = 0.5, angle = 0))
p <- p + geom_boxplot(data = proportion_plot, aes(x = Group, y = Value, fill = Group), outlier.shape = NA)
p <- p + geom_point(data = proportion_plot, aes(x = Group, y = Value, fill = Group), size=0.8, position = position_jitter(width = 0.25, height = 0)) 
p <- p + facet_wrap(~Celltype, ncol = 5, scales = "free_y") +
  theme(
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines")
  )
print(p)
ggsave('cell_type_composition_primary_type_boxplot.pdf', p, width = 12, height = 14)

##################### 箱线图 proportion figure  by    samplename   ############################
library(ggplot2)

# 绘制箱线图和小提琴图
p <- ggplot(proportion_plot, aes(x = Group, y = Value, fill = Group)) +
  # 添加箱线图
  geom_boxplot(outlier.shape = NA, width = 0.7, alpha = 0.7) +  # 隐藏异常值点，设置宽度和透明度
  # 添加散点图
  geom_point(
    size = 0.8, 
    position = position_jitter(width = 0.2, height = 0),  # 调整点的分布
    alpha = 0.6,  # 设置点的透明度
    color = "black"  # 设置点的颜色
  ) +
  # 分面显示
  facet_wrap(
    ~Celltype, 
    ncol = 5,  # 每行显示 5 个分面
    scales = "free_y"  # y 轴独立缩放
  ) +
  # 设置主题
  theme_minimal(base_size = 12) +  # 使用简洁主题，设置基础字体大小
  theme(
    panel.background = element_rect(fill = "white"),  # 设置背景为白色
    axis.line = element_line(color = "black", size = 0.15),  # 设置坐标轴线
    axis.text.x = element_text(
      size = 7, 
      color = "black", 
      vjust = 0.5, 
      hjust = 0.5, 
      angle = 90  # 设置 x 轴文字角度
    ),
    axis.text.y = element_text(
      size = 7, 
      color = "black", 
      vjust = 0.5, 
      hjust = 0.5
    ),
    strip.placement = "outside",  # 分面标题放置在外部
    strip.text = element_text(size = 10, face = "bold"),  # 设置分面标题样式
    panel.spacing = unit(1.5, "lines"),  # 调整分面图间距
    legend.position = "right",  # 图例放置在右侧
    legend.title = element_text(size = 10, face = "bold"),  # 图例标题样式
    legend.text = element_text(size = 8)  # 图例文字样式
  ) +
  # 设置坐标轴标签
  labs(
    x = "Group", 
    y = "Percentage (%)", 
    fill = "Group"
  )

# 打印图形
print(p)

# 保存图像
ggsave(
  filename = "cell_type_composition_primary_type_boxplot.pdf", 
  plot = p, 
  width = 12, 
  height = 14, 
  dpi = 300
)

##########横向堆积图 by group############
legend_order <- desired_primary_order

# 转换因子顺序
proportion_plot$Group <- factor(
  proportion_plot$Group,
  levels = group_order  # 强制分组顺序
)

proportion_plot$Celltype <- factor(
  proportion_plot$Celltype,
  levels = legend_order  # 强制图例顺序
)

# 确保颜色顺序与 legend_order 一致
my_cols_primary <- my_cols_primary[legend_order]

# 数据预处理（关键步骤）
group_proportion <- proportion_plot %>%
  group_by(Group, Celltype) %>%
  summarise(Value = sum(Value), .groups = "drop") %>%  # 按组汇总比例
  group_by(Group) %>%
  mutate(Total = sum(Value),
         Percent = Value/Total*100) %>%  # 标准化到100%
  ungroup()

# 可视化代码
p <- ggplot(group_proportion, 
            aes(x = Group, y = Percent, fill = Celltype)) +
  geom_col(position = "stack", width = 0.7, 
           color = "white", linewidth = 0.2, key_glyph = "rect") +
  coord_flip() +
  scale_fill_manual(values = my_cols_primary, 
                    breaks = legend_order,
                    guide = guide_legend(
                      reverse = TRUE,
                      keywidth = 0.5,  # 关键调整：压缩图例色块宽度
                      keyheight = 0.8   # 保持色块高度避免变形
                    )) +
  scale_y_continuous(expand = c(0, 0), 
                     labels = scales::label_number(suffix = "%")) +
  labs(
    title = "Cell Type Composition Across Experimental Groups",
    x = NULL,
    y = NULL,
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(color = "black", size = 10, family = "Arial"),
    axis.text.y = element_text(color = "black", size = 12, face = "bold", margin = margin(r = 10)),  # 减少右侧边距
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),  # 进一步缩小图例键
    legend.spacing.y = unit(0.05, "cm"),  # 最小化图例间距
    legend.title = element_text(face = "bold", size = 9),  # 缩小标题
    legend.text = element_text(size = 8),  # 缩小文本
    legend.margin = margin(0, 0, 0, -5),  # 压缩图例右侧边距
    plot.margin = margin(1, 1.5, 1, 1, "cm"),  # 减少右侧整体边距
    axis.line.x = element_line(color = "black", linewidth = 0.5)
  )
# 输出设置
ggsave(
  "横向cell_type_composition_primary_type_stacked_bar_chart_by group.pdf", 
  plot = p,
  width = 13,  # 显著增加宽度
  height = 4,   # 按比例调整高度
  device = cairo_pdf,
  dpi = 300
)

#####横向堆积图 by samplename########
p <- ggplot(proportion_plot, aes(x = reorder(Name, Value), y = Value, fill = Celltype)) +
  geom_col(position = "stack", width = 0.85, color = "white", linewidth = 0.1) +  # 优化条形外观
  coord_flip() +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +  # 垂直分面
  scale_fill_manual(values = my_cols_primary, breaks = legend_order) +
  scale_y_continuous(expand = c(0, 0)) +  # 去除y轴空白
  labs(
    title = "Cell Type Composition by Group",
    x = NULL,
    y = "Percentage (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9, color = "black", margin = margin(r = 5)),
    axis.title.x = element_text(size = 12, face = "plain", margin = margin(t = 8)),
    strip.text.y = element_text(angle = 0, size = 12, face = "bold", hjust = 0),  # 分面标签左对齐
    strip.placement = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
    panel.spacing.y = unit(1.5, "lines"),  # 增加分面间距
    legend.position = "right",
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin = margin(1, 1.5, 1, 1, "cm")
  ) +
  guides(fill = guide_legend(byrow = TRUE, title.position = "top"))

# 动态计算图片高度（基于分组和样本数量）
group_counts <- proportion_plot %>% 
  distinct(Group, Name) %>% 
  count(Group) 

base_height <- 6  # 基础高度
height_per_sample <- 0.3  # 每个样本的高度系数
total_height <- base_height + sum(group_counts$n) * height_per_sample

# 保存图像
ggsave(
  filename = "横向cell_type_composition_primary_type_stacked_bar_chart_by samplename.pdf",
  plot = p,
  width = 12,  # 固定宽度
  height = min(total_height, 50),  # 安全高度限制
  dpi = 300,
  limitsize = FALSE
)
