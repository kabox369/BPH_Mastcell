library(Seurat)
library(rjson)
library(dplyr)
source("/home/cao/公共的/Tang/configs/color.R")
source("~/Single Cell/configs/selfFunction.R")

meta_list <- list()
meta_colname_list <- list()
primary_type_list <- c()

file_list <- list.files("~/Single Cell/BPH_LS/secondary annotation/secondary annotation RDS")
for (i in file_list) {
  cat(paste("Read meta: <", i, "> is running...\n"))
  # define constant
  meta_file_path <- file.path("~/Single Cell/BPH_LS/secondary annotation/secondary annotation RDS", i)
  
  # read meta and remove columns name: cluster.X.X.X
  obj <- readRDS(meta_file_path)
  meta_cache <- obj@meta.data
  cluster_clean_list <- colnames(meta_cache)[grepl("clusters", colnames(meta_cache))]
  meta_cache[, cluster_clean_list] <- NULL
  
  meta_list[[i]] <- meta_cache
  meta_colname_list[[i]] <- colnames(meta_cache)
}

new_meta <- data.frame()
keep_colnames <- Reduce(intersect, meta_colname_list)
for (i in file_list) {
  new_meta <- rbind(new_meta, meta_list[[i]][, keep_colnames])
}

unique(new_meta$secondary_type)

seurat_object <- readRDS("~/Single Cell/BPH_LS/BPH_LS.major.annotation.pc40umap30res0.8.rds")
seurat_object$old_primary_type <- seurat_object$primary_type

# new_meta$secondary_type <- stringi::stri_replace_all_fixed(new_meta$secondary_type, "+", "")
# new_meta$secondary_type <- stringi::stri_replace_all_fixed(new_meta$secondary_type, " ", "")
# new_meta$secondary_type <- stringi::stri_replace_all_fixed(new_meta$secondary_type, "-", "")
# new_meta$secondary_type <- stringi::stri_replace_all_fixed(new_meta$secondary_type, "_", "")
# new_meta$secondary_type <- stringi::stri_replace_all_fixed(new_meta$secondary_type, "?", "1")

mapInfo <- rjson::fromJSON(file = '~/Single Cell/BPH_LS/nameing json/map.json')
json_filepath <- '~/Single Cell/BPH_LS/nameing json/map.json'
json_content <- readChar(json_filepath, file.info(json_filepath)$size)
mapInfo <- fromJSON(json_content)


setdiff(names(mapInfo), unique(new_meta$secondary_type))
setdiff(unique(new_meta$secondary_type), names(mapInfo))


pass_cell <- setdiff(rownames(seurat_object@meta.data), rownames(new_meta))
new_meta_passcell <- seurat_object@meta.data[pass_cell, intersect(colnames(seurat_object@meta.data), colnames(new_meta))]
new_meta_passcell$secondary_type <- "PASS"

convert_all_columns_to_character <- function(data) {
  for (col_name in names(data)) {
    data[[col_name]] <- as.character(data[[col_name]])
  }
  return(data)
}

data1 <- new_meta[intersect(rownames(seurat_object@meta.data), rownames(new_meta)), ]
data2 <- new_meta_passcell[, colnames(new_meta)]
data1 <- convert_all_columns_to_character(data1)
data2 <- convert_all_columns_to_character(data2)
new_meta_add <- bind_rows(data1, data2)

length(seurat_object@meta.data$orig.ident) == length(new_meta_add$orig.ident)
seurat_object <- AddMetaData(seurat_object, new_meta_add)
seurat_object <- label_clusters(seurat_object, mapInfo)

table(seurat_object@meta.data$primary_type)
table(seurat_object@meta.data$secondary_type)
DimPlot(seurat_object, reduction = "harmony.major.pcs40.umap.dims30", group.by = "primary_type", label = F, 
        cols = my_cols_primary, alpha = 0.2, pt.size = 0.2 ,raster.dpi = c(1080,1080) ,raster = F)
DimPlot(seurat_object, reduction = "harmony.major.pcs40.umap.dims30", group.by = "secondary_type", label = T,
        cols = my_cols_secondary, alpha = 0.2, pt.size = 0.2 ,raster.dpi = c(1080,1080) ,raster = F)

seurat_object <- subset(seurat_object,subset = primary_type != "PASS")
saveRDS(seurat_object, '~/Single Cell/BPH_LS/BPH_LS_mappedLabelled.rds',compress = F)

