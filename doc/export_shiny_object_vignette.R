## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(cellcuratoR)
library(Seurat)
library(tidyverse)

if(!requireNamespace("SeuratData", quietly = TRUE)){
  stop("Package \"SeuratData\" needed for this vignette to work. Please install it.",
      call. = FALSE)
}

library(SeuratData)

InstallData("ifnb") 
data("ifnb")

## ------------------------------------------------------------------------
head(ifnb@meta.data)

## ------------------------------------------------------------------------
ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# We map each cluster to a putative cell type
seurat_celltypes <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", 
                      "CD16 Mono", "B", "CD8 T", "T activated", 
                      "NK", "DC", "B Activated", "Mk",  "pDC", "Eryth")
celltype <- plyr::mapvalues(x = immune.combined@meta.data$seurat_clusters, 
                                from = seq(from = 0, to = 12), 
                                to = seurat_celltypes)
immune.combined@meta.data$celltype <- celltype

# We add a treatment column to the seurat object so that differential expression can 
# be performed between unstimulated and stimulated samples. Treatment must be a binary factor. 
treatment <- as.factor(immune.combined@meta.data$orig.ident)
immune.combined@meta.data <- data.frame(immune.combined@meta.data, treatment)

## ------------------------------------------------------------------------
## consider changing filepath to ~/Desktop or other directory that is easier to navigate to
my_filepath <- file.path(find.package("cellcuratoR"))
dir.create(file.path(my_filepath, "my_cellcuratoR_objects/"))


## ----message = FALSE, warning = FALSE------------------------------------
export_shiny_object(seurat_object = immune.combined, 
                    final_cluster_column_name = "seurat_clusters",
                    library_id_column_name = "orig.ident",
                    classification_column_name = "celltype",
                    create_dendrogram = TRUE,
                    custom_cluster_colors = FALSE,
                    custom_library_colors = FALSE,
                    additional_metadata_cols = c("treatment"), # allows for dge between treatment groups
                    export_data_path = file.path(my_filepath, "my_cellcuratoR_objects/"))
## consider changing filepath to ~/Desktop or other directory that is easier to navigate to

## ------------------------------------------------------------------------
print(paste0("exporting data to: ", my_filepath, "my_cellcuratoR_objects/"))
# Launch the app with the following command (commented out for R-markdown)
# cellcuratoR::launchApp()

# 1. Navigate to the our newly created directory printed above (or change to a more accesible directory for you)
# 2. Select the "infb_shiny" dataset from the "which dataset should be loaded?" dropdown
# 3. Interactively explore data!


