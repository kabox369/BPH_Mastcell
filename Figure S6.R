#Figure S6A-B
go_enrichment("~/Single Cell/BPH_LS/secondary cell type figure/t cell",gene_list = T_DEG$gene,species = "human",ont = "BP", padj = "BH")
go_enrichment("~/Single Cell/BPH_LS/secondary cell type figure/b cell",gene_list = B_DEG$gene,species = "human",ont = "BP", padj = "BH")

write.csv(GO_enrich_res@result,"B cell GO pathway.csv")
write.csv(KEGG_enrich_res@result,"B cell KEGG pathway.csv")

dotplot(GO_enrich_res, showCategory = B_cell_GO_selected$Description, x = "GeneRatio",color = "p.adjust",title = NULL)
dotplot(KEGG_enrich_res, showCategory = T_cell_KEGG$Description, x = "GeneRatio",color = "p.adjust",title = NULL) 

#Figure S6C
object.list <- list(Large = BPH_Large_EndoMT_Cellchat, 
                    Small = BPH_Small_EndoMT_Cellchat )
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked =T, do.stat = TRUE,x.rotation = 0,
        comparison = c(1, 2),color.use = c("#2CA02C","#FF7F0E"),)