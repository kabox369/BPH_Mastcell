#Figure S11 A-B 
par(mfrow=c(1,1))
netVisual_aggregate(BPH_EndoMT_Cellchat, signaling = "TGFb", layout = "circle", signaling.name = 'TGF-β Signaling Pathways', 
                    targets.use = cell_intrest,top = 0.2,vertex.weight.max = 1)

par(mfrow=c(1,1))
netVisual_aggregate(BPH_EndoMT_Cellchat, signaling = "VEGF", layout = "circle", signaling.name = 'VEGF Signaling Pathways', 
                    targets.use = cell_intrest,top = 0.2,vertex.weight.max = 1)

#Figure S11C
Idents(object_mast)="secondary_type"
custom_order <- c("Mast cell_CDKN1A","Mast cell_JAK2","Mast cell_CTSG")
gene_list <- c("CDKN1A","NFKB1","DUSP2","VEGFA",
               "JAK2","RUNX1","DAPK1",
               "CTSG","COMMD6"
)
plot_df <- FetchData(object_mast, vars = c(gene_list, "secondary_type"), layer = "data") %>%
  mutate(
    secondary_type = factor(secondary_type, levels = custom_order)  
  ) %>%
  group_by(secondary_type) %>%
  summarise_all(function(x) mean(x = expm1(x = x))) %>%
  filter(!is.na(secondary_type)) %>%  
  arrange(secondary_type) %>%  
  column_to_rownames("secondary_type") %>%
  log2()
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#0F7B9F", "white", "#D83215")
)
plot_df_t <- plot_df %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>% 
  mutate(Gene = factor(Gene, levels = gene_list)) %>%  
  column_to_rownames("Gene") %>% 
  as.matrix()

hm <- ComplexHeatmap::Heatmap(
  plot_df_t,  
  col = col_fun,
  name = "Gene expression",
  
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 12, fontface = "italic", col = "black"),
  column_names_gp = gpar(fontsize = 12),
  
  rect_gp = gpar(col = "gray90", lwd = 0.5),  
  border_gp = gpar(col = "black", lwd = 1.5), 
  
  heatmap_legend_param = list(
    legend_direction = "vertical",
    title_position = "topcenter",
    legend_height = unit(2, "cm"),
    legend_width = unit(1, "cm"),
    gap = unit(0.5, "cm"),
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  right_annotation = rowAnnotation(
    spacer = anno_empty(
      width = unit(2, "cm"),  
      border = FALSE
    )
  ),
  width = ncol(plot_df_t)*unit(12, "mm"),  
  height = nrow(plot_df_t)*unit(10, "mm")  
)
pdf("mast_vertical_heatmap.pdf", width = 10, height = 16, useDingbats = FALSE)

draw(
  hm,
  padding = unit(c(2, 4, 2, 2), "cm"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  gap = unit(2, "cm"),
  adjust_annotation_extension = TRUE
)

dev.off()

#Figure S11 E
#aftre mast cell monocle3 analyse
suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(Matrix)
  library(scales)
})

genes <- c("TPSAB1", "FCER1A", "VEGFA", "TGFB1")

# cds must already contain UMAP in reducedDims(cds)$UMAP
stopifnot("UMAP" %in% names(reducedDims(cds)))

# keep genes that exist
genes_use <- intersect(genes, rowData(cds)$gene_short_name %||% rownames(cds))
if (length(genes_use) == 0) stop("None of the requested genes found in cds.")

# map gene symbols to rownames if gene_short_name exists
rn <- rownames(cds)
gshort <- rowData(cds)$gene_short_name
if (!is.null(gshort)) {
  idx <- match(genes_use, gshort)
  idx <- idx[!is.na(idx)]
  rn_use <- rn[idx]
  genes_use <- gshort[idx]
} else {
  rn_use <- genes_use
}

um <- as.data.frame(reducedDims(cds)$UMAP)
colnames(um) <- c("UMAP_1", "UMAP_2")
um$cell <- rownames(um)

# expression matrix (counts); change to "logcounts" if you have it
expr_mat <- assays(cds)[["counts"]]
expr_mat <- expr_mat[rn_use, colnames(cds), drop = FALSE]

# build long dataframe
df_list <- vector("list", length(rn_use))
for (i in seq_along(rn_use)) {
  g <- genes_use[i]
  v <- as.numeric(expr_mat[i, ])
  df_list[[i]] <- data.frame(
    cell = colnames(cds),
    gene = g,
    expr = v,
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, df_list)
df <- merge(df, um, by = "cell", all.x = TRUE)

# log10(Expression); color only for expr > 0
df$expr_log10 <- log10(df$expr + 1)
df$expr_log10[df$expr == 0] <- NA

p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey80", size = 1.1) +
  geom_point(
    data = df[!is.na(df$expr_log10), ],
    aes(color = expr_log10),
    size = 1.1
  ) +
  facet_wrap(~ gene, ncol = 2) +
  scale_color_viridis_c(
    name = "log10(Expression)",
    limits = c(0, 2),
    oob = squish
  ) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.position = "right"
  )

p

#Figure S11 F
Idents(object_mast)="group"
custom_order <- c("Large","Small")
gene_list <- c("VEGFA","NAMPT","TGFB1","TNFSF10","TULP2")
plot_df <- FetchData(object_mast, vars = c(gene_list, "group"), layer = "data") %>%
  mutate(
    secondary_type = factor(secondary_type, levels = custom_order)  
  ) %>%
  group_by(secondary_type) %>%
  summarise_all(function(x) mean(x = expm1(x = x))) %>%
  filter(!is.na(secondary_type)) %>%  
  arrange(secondary_type) %>%  
  column_to_rownames("secondary_type") %>%
  log2()
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#0F7B9F", "white", "#D83215")
)
plot_df_t <- plot_df %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>% 
  mutate(Gene = factor(Gene, levels = gene_list)) %>%  
  column_to_rownames("Gene") %>% 
  as.matrix()

hm <- ComplexHeatmap::Heatmap(
  plot_df_t,  
  col = col_fun,
  name = "Gene expression",
  
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 12, fontface = "italic", col = "black"),
  column_names_gp = gpar(fontsize = 12),
  
  rect_gp = gpar(col = "gray90", lwd = 0.5),  
  border_gp = gpar(col = "black", lwd = 1.5), 
  
  heatmap_legend_param = list(
    legend_direction = "vertical",
    title_position = "topcenter",
    legend_height = unit(2, "cm"),
    legend_width = unit(1, "cm"),
    gap = unit(0.5, "cm"),
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  right_annotation = rowAnnotation(
    spacer = anno_empty(
      width = unit(2, "cm"),  
      border = FALSE
    )
  ),
  width = ncol(plot_df_t)*unit(12, "mm"),  
  height = nrow(plot_df_t)*unit(10, "mm")  
)
pdf("mast_vertical_heatmap_by_group.pdf", width = 10, height = 16, useDingbats = FALSE)

draw(
  hm,
  padding = unit(c(2, 4, 2, 2), "cm"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  gap = unit(2, "cm"),
  adjust_annotation_extension = TRUE
)

dev.off()

#########Figure S11D 
p4 <- DimPlot(
  object_mast , 
  reduction = "UMAP", 
  group.by = "secondary_type", 
  label = FALSE, 
  cols = my_cols_secondary, 
  alpha = 1, 
  pt.size = 1, 
  raster = F,
  split.by = "group"
) + ggtitle(NULL)

ggsave("Mast cell secondary annotation figure split by group.pdf", p4, width = 12, height = 6)
