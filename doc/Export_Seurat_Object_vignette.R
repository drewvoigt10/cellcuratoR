## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----cars----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(SCOTA)

head(pbmc_small@meta.data)

## ------------------------------------------------------------------------
celltype <- plyr::mapvalues(x = pbmc_small@meta.data$RNA_snn_res.1, 
                            from = seq(from = 0, to = 2), 
                            to = c("celltype_A", "celltype_B", "celltype_C"))
pbmc_small@meta.data <- data.frame(pbmc_small@meta.data, celltype)

## ----message = FALSE, warning = FALSE------------------------------------
export_shiny_object(seurat_object = pbmc_small, 
                    final_cluster_column_name = "RNA_snn_res.1",
                    library_id_column_name = "orig.ident",
                    cell_classification_column_name = "celltype",
                    create_dendrogram = TRUE,
                    custom_cluster_colors = FALSE,
                    custom_library_colors = FALSE,
                    export_data_path = "~/Desktop/pbmc_small_shiny/")

