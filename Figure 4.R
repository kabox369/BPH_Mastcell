#Figure 4A
library(CellChat)
library(ggplot2)
library(patchwork)

cc = BPH_LS_mappedLabelled_EndoMT_Cellchat
p_net <- netVisual_aggregate(
  object = cc,
  signaling = "KIT",
  layout = "circle",
  edge.weight.max = NULL,  # let it scale automatically
  vertex.weight = NULL,    # or: as.numeric(table(cc@idents))
  title.name = "KIT signaling pathway network"
)

# Option A: let CellChat choose pathway-related genes
p_vln <- plotGeneExpression(
  object = cc,
  signaling = "KIT",     # uses pathway genes
  type = "violin",
  show.legend = FALSE
)

# Option B: force exactly the two genes shown in your figure
# p_vln <- plotGeneExpression(
#   object = cc,
#   gene.use = c("KITLG", "KIT"),
#   type = "violin",
#   show.legend = FALSE
# )

p <- p_net / p_vln + plot_layout(heights = c(2.2, 1))
p
