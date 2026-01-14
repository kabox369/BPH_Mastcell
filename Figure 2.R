#Figure 2A
suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(Matrix)
  library(scales)
  library(patchwork)
})

reduction <- "harmony.pcs40.umap.dims30"

cds <- as.cell_data_set(object_fibro)

umap_mat <- Embeddings(object_fibro, reduction = reduction)
umap_mat <- umap_mat[colnames(cds), , drop = FALSE]
reducedDims(cds)$UMAP <- umap_mat

colData(cds)$seurat_clusters <- Idents(object_fibro)[colnames(cds)] |> as.character()

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

p_traj <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "seurat_clusters",
  label_groups_by_cluster = TRUE,
  label_branch_points = TRUE,
  label_leaves = TRUE
)
p_traj
if (!"UMAP" %in% names(reducedDims(cds))) stop("UMAP not found in reducedDims(cds)$UMAP")
if (is.null(principal_graph(cds)[["UMAP"]])) stop("No principal graph found. Run learn_graph(cds) first.")
if (!"pseudotime" %in% names(colData(cds))) stop("No pseudotime found. Run order_cells(cds, root_cells=...) first.")

plot_gene_umap_log10 <- function(cds, gene, assay_name = "counts",
                                 lim = c(0, 2),
                                 pt_size = 0.6) {
  um <- as.data.frame(reducedDims(cds)$UMAP)
  colnames(um) <- c("UMAP_1", "UMAP_2")
  um$cell <- rownames(um)
  
  gshort <- rowData(cds)$gene_short_name
  rn <- rownames(cds)
  
  if (!is.null(gshort) && gene %in% gshort) {
    ridx <- which(gshort == gene)[1]
    rnm <- rn[ridx]
  } else if (gene %in% rn) {
    rnm <- gene
  } else {
    stop("Gene not found in cds: ", gene)
  }
  
  if (!assay_name %in% names(assays(cds))) {
    stop("Assay not found in cds assays: ", assay_name,
         ". Available: ", paste(names(assays(cds)), collapse = ", "))
  }
  
  expr_mat <- assays(cds)[[assay_name]]
  v <- as.numeric(expr_mat[rnm, colnames(cds)])
  df <- data.frame(cell = colnames(cds), expr = v, stringsAsFactors = FALSE)
  df <- merge(df, um, by = "cell", all.x = TRUE)
  
  df$expr_log10 <- log10(df$expr + 1)
  df$expr_log10[df$expr == 0] <- NA
  
  ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "grey80", size = pt_size) +
    geom_point(
      data = df[!is.na(df$expr_log10), ],
      aes(color = expr_log10),
      size = pt_size
    ) +
    scale_color_viridis_c(
      name = "log10(Expression)",
      limits = lim,
      oob = squish
    ) +
    labs(title = gene, x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5),
      legend.position = "right"
    )
}

p1 <- plot_gene_umap_log10(cds, "PECAM1", assay_name = "counts", lim = c(0, 1.2), pt_size = 0.7)
p2 <- plot_gene_umap_log10(cds, "ACTA2",  assay_name = "counts", lim = c(0, 2.0), pt_size = 0.7)

p3 <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.8,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  label_cell_groups = FALSE
) +
  theme_classic(base_size = 11)

(p1 | p2 | p3)

#Figure 2B
DotPlot(object_fibro_EndMT,features = c('VWF','PECAM1','ACTA2',
                                   'MCAM','CD34','KDR',
                                   'PTPRC','CD14'), min.cutoff = 0.5, cols = c("gray", "coral2"))
#Figure 2H
library(dplyr)
library(ggplot2)
library(forcats)
library(scales)
cc = BPH_Large_EndoMT_Cellchat 
immune_cells <- c("Myeloid cell","T cell","B cell","NK","Mast cell","Dendritic cell")
key_lig_chemo  <- c("CXCL12","CCL2","CCL5","MIF","IL6","CSF1","EDN1")
key_lig_growth <- c("FGF2","FGF7","HGF","PDGFA","PDGFB","IGF1","EGF","ANGPTL4","PTN","MDK")
key_lig_ecm    <- c("COL1A1","COL1A2","COL4A1","COL4A2","FN1","COL6A2")

# -------- Step 1: Raw data -------- #
df_raw <- subsetCommunication(cc) %>%
  filter(
    source == "Stromal cell_PECAM1+ACTA2+",
    target %in% c("Fibroblast", immune_cells),
    ligand %in% c(key_lig_chemo, key_lig_growth, key_lig_ecm)
  ) %>%
  mutate(
    target = ifelse(target %in% immune_cells, "Immune", "Fibroblast"),
    class = case_when(
      ligand %in% key_lig_chemo  ~ "Chemo",
      ligand %in% key_lig_growth ~ "Growth",
      ligand %in% key_lig_ecm    ~ "ECM"
    )
  ) %>%
  mutate(target = factor(target, levels = c("Fibroblast","Immune")))

# -------- Step 2: Summaries -------- #
bubble <- df_raw %>%
  group_by(class, ligand, target) %>%
  summarise(score = sum(prob), .groups = "drop") %>%
  group_by(class, target) %>%
  mutate(
    score_rel = score / sum(score),
    label_pct = percent(score_rel, accuracy = 1)
  ) %>% 
  ungroup()

# -------- Step 3: ordering -------- #
bubble$ligand <- fct_reorder(bubble$ligand, bubble$score)

# log-scale point size
bubble$size_scaled <- log10(bubble$score + 1e-6) %>% rescale(to = c(1, 8))

# -------- Step 4: Plot with BOTH-SIDE percentages -------- #
p <- ggplot(bubble, aes(x = target, y = ligand)) +
  geom_point(aes(size = size_scaled, fill = class),
             shape = 21, color = "black", stroke = 0.25) +
  
  geom_text(
    data = subset(bubble, target == "Fibroblast"),
    aes(x = as.numeric(target) + 0.25, label = label_pct),
    size = 3.2, hjust = 0
  ) +
  
  geom_text(
    data = subset(bubble, target == "Immune"),
    aes(x = as.numeric(target) + 0.25, label = label_pct),
    size = 3.2, hjust = 0
  ) +
  
  scale_size_identity() +
  scale_fill_manual(values = c(
    Chemo="#56B4E9",
    ECM="#9B59B6",
    Growth="#E69F00"
  )) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(clip="off") +
  labs(
    title = "Key soluble + ECM ligands from PECAM1⁺ ACTA2⁺ stromal cell",
    subtitle = "Dot size = absolute strength (log-scaled) | Side labels = % within each class",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size=12, face="bold"),
    legend.position = "none",
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=10),
    plot.title = element_text(size=14, face="bold", hjust=0.51),
    plot.subtitle = element_text(size=12, hjust=0.51),
    plot.margin = margin(8, 80, 8, 80)
  )

p