#' Opposite of in
#'
#' \code{%!in%} returns a vector of the positions of the first vector that are
#' not in the second vector. Adapted from stackoverflow Sacha Epskamp
#' https://stackoverflow.com/questions/5831794/opposite-of-in
#'
#' @param x The first vector.
#'
#' @param y The second vector.
#'
#' @example
#' \dontrun{
#' c(2,3,4,5) %!in% c(3,4)
#'}
"%!in%" <- function(x, y) ! ("%in%"(x, y))



#' Emulation of gg colors
#'
#' \code{gg_color_hue} Returns a vector containing n unique hex colors
#' Adapted from stackoverflow John Colby
#' https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

#' @param n The number of colors to return.
#'
#' @example
#' \dontrun{
#' gg_color_hue(10)
#'}


gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Export seurat-based shiny object
#'
#' \code{export_shiny_object} returns a stripped down seurat object, a
#' dendrogram, colors, and a list of mapped genes for input into the cellcuratoR
#' shiny application. This function aims to help researchers shiny-ize their
#' own datasets for local interactive analysis.
#'
#' @param seurat_objectIinput S4 Seurat object.
#'
#' @param final_cluster_column_name The column name within the meta.data
#' (character string) that corresponds to the final cluster label.
#'
#' @param library_id_column_name The column name within the meta.data
#' (character string) that corresponds to the library that each cell belongs to.
#' If CCA was performed to aggregate samples from multiple libraries, this
#' is usually stored in the column "orig.ident".
#'
#' @param classification_column_name The column name within the meta.data
#' (character string) that corresponds to the cell type that each
#' cell belongs to. If the user has not manually mapped each cell to a cell
#' type, simply provide the column used for the library_id_column_name
#' argument. Cell types must perfectly correspond to the argument
#' final_cluster_column_name (eg, all cells in cluster 1 must have the
#' same cell type classification)
#'
#' @param create_dendrogram TRUE/FALSE argument if a dendrogram should be
#' created. The dendrogram is used by the shiny interface to organize the
#' appearance of clusters in violin plots. User should select TRUE unless
#' they have built a dendrogram with custom settings. Defaults to TRUE.
#'
#' @param replacement_dendrogram If the user has already created a
#' dendrogram of cell clusters, the dendrogram can be supplied as a
#' .RData object. Defaults to NULL.
#'
#' @param custom_cluster_colors If the user would like output plots to have
#' custom colors for the identity final_cluster_column_name, a vector
#' consisting  of colors (with length = n_clusters) can be supplied.
#' Defaults to FALSE.
#'
#' @param custom_library_colors If the user would like output plots to have
#' custom colors for the identity library_id_column_name, a vector consisting
#' of colors (with length = n_libraries) can be supplied. Otherwise, specify
#' FALSE. Defaults to FALSE.
#'
#'@param additional_metadata_cols Additional columns in the meta.data can be
#' supplied that are available for differential expression in the shiny
#' interface. Columns must be binary factors. Deafults to NULL.
#'
#' @param export_data_path A character string of the directory path where data
#' should be exported. Of note, all exported data objects for shiny should be
#' stored in the same directory in order to be accessible at the same time within
#' the shiny application.
#'
#' @examples
#' \dontrun{
#' celltype <- plyr::mapvalues(x = pbmc_small@meta.data$RNA_snn_res.1,from = seq(from = 0, to = 2), to = c("celltype_A", "celltype_B", "celltype_C"))
#' pbmc_small@meta.data <- data.frame(pbmc_small@meta.data, celltype)
#' pbmc_small@meta.data$celltype <- celltype
#' export_shiny_object(seurat_object = pbmc_small,
#' final_cluster_column_name = "RNA_snn_res.1",
#' library_id_column_name = "orig.ident",
#' classification_column_name = "celltype",
#' create_dendrogram = TRUE,
#' custom_cluster_colors = FALSE,
#' custom_library_colors = FALSE,
#' export_data_path = "~/Desktop/pbmc_small_shiny/")
#' }
#'
#' @export
export_shiny_object <- function(seurat_object,
                                final_cluster_column_name,
                                library_id_column_name,
                                classification_column_name,
                                create_dendrogram = TRUE,
                                replacement_dendrogram = NULL,
                                custom_cluster_colors = FALSE,
                                custom_library_colors = FALSE,
                                additional_metadata_cols = NULL,
                                export_data_path) {
  DefaultAssay(seurat_object) <- "RNA" # default assay of
                                       # object must be RNA
  # create export directory
  dir.create(export_data_path, recursive = TRUE)

  #############################################
  ##### strip off data from seurat object #####
  #############################################

  my_reductions <- seurat_object@reductions
  my_metadata <- seurat_object@meta.data
  my_version <- seurat_object@version
  my_assays <- seurat_object@assays
  my_active_assay <- seurat_object@active.assay
  my_active_ident <- seurat_object@active.ident

  if (final_cluster_column_name %!in% colnames(my_metadata) |
     library_id_column_name %!in% colnames(my_metadata) |
     classification_column_name %!in% colnames(my_metadata)) {

    name_vector <- c(final_cluster_column_name,
                    library_id_column_name,
                    classification_column_name)

    stop(paste(name_vector[
      which(
        name_vector %!in% colnames(my_metadata)
        )
      ],
      "not a column name in the meta.data",
      sep = " "))

  } else {
    if (is.null(additional_metadata_cols)) {
      my_final_metadata <- my_metadata[,
                                       c(final_cluster_column_name,
                                         library_id_column_name,
                                         classification_column_name)
                                       ]
      colnames(my_final_metadata) <- c("final_cluster_labels",
                                       "libraryID",
                                       "celltype")
    } else {
      my_final_metadata <- my_metadata[,
                                       c(final_cluster_column_name,
                                         library_id_column_name,
                                         classification_column_name,
                                         additional_metadata_cols)
                                       ]
      colnames(my_final_metadata)[1:3] <- c("final_cluster_labels",
                                            "libraryID",
                                            "celltype")
    }

  }

  #############################################
  #####  create colors for the umap plot  #####
  #############################################

  # Colors for the tSNE plot are based on the number of unique elements
  # in the final_cluster_labels column
  n_clusters <- length(
    summary(
      as.factor(
        my_final_metadata$final_cluster_labels
        )
      )
    )
  # if a vector is supplied to custom_cluster_colors
  if (length(custom_cluster_colors) > 1) {
    if (length(custom_cluster_colors) == n_clusters) {
      final_colors <- custom_cluster_colors
    } else {
      stop("length of supplied custom colors not the same as the
           number of clusters in the column final_cluster_column_name")
    }
  # if the user did not supply a custom color vector, we create a color
  # vector for them with gg_color_hue
  } else {
    n_colors <- ceiling(n_clusters / 2) # number of total colors to create
    dark <- gg_color_hue(n_colors)     # dark colors
    pastel <- colorspace::qualitative_hcl(n_colors, "Pastel 1") # light colors

    final_colors <- c()
    for (i in 1:n_colors) {
      dark_col <- dark[i]
      pastel_col <- pastel[i]
      final_colors <- c(final_colors, dark_col, pastel_col)
    }

    final_colors <- final_colors[1:n_clusters]
  }

  # Colors for the libraryID are based on the number of unique elements in the
  # libraryID column
  n_libraries <- length(summary(as.factor(my_final_metadata$libraryID)))
  # if a vector is supplied to custom_cluster_colors
  if (length(custom_library_colors) > 1) {
    if (length(custom_library_colors) == n_libraries) {
      final_library_colors <- custom_library_colors
    } else {
      stop("length of supplied custom colors not the same as the number of
           clusters in the column library_id_column_name")
    }
  # if the user did not supply a custom library color vector, we create a
  # library color vector for them with gg_color_hue
  } else {
    n_colors <- ceiling(n_libraries / 2) # number of total colors to create
    dark <- gg_color_hue(n_colors)     # dark colors
    pastel <- colorspace::qualitative_hcl(n_colors, "Pastel 1") # light colors

    final_library_colors <- c()
    for (i in 1:n_colors) {
      dark_col <- dark[i]
      pastel_col <- pastel[i]
      final_library_colors <- c(final_library_colors, dark_col, pastel_col)
    }

    final_library_colors <- final_library_colors[1:n_libraries]
  }


  #############################################
  ##### dendrogram export (if requested)  #####
  #############################################

  # create a dendrogram based on clusters (final_cluster_labels) that
  # express similar genes we will use the original seurat object as it
  # still has expression data for all genes we will identify the 100
  # most enriched genes in each cluster before hierarchical clustering
  if (create_dendrogram == TRUE) {
    Idents(seurat_object) <- final_cluster_column_name

    identity_vector <- names(
      summary(
        as.factor(
          seurat_object@active.ident
          )
        )
      )

    all_enriched_genes <- c() # stores top 100 enriched genes in each cluster
    for (i in 1:n_clusters) {
      de_genes <- FindMarkers(seurat_object,
                              ident.1 = identity_vector[i],
                              min.pct = 0.5,
                              logfc.threshold = 0.5)

      enriched_genes <- de_genes %>%
        rownames_to_column("gene") %>%
        filter(avg_logFC > 0) %>%
        arrange(desc(avg_logFC))

      if (nrow(enriched_genes) > 100) {
        all_enriched_genes <- c(all_enriched_genes, enriched_genes[1:100, 1])
      } else {
        all_enriched_genes <- c(all_enriched_genes, enriched_genes[, 1])
      }
    }
    final_enriched_genes <- unique(all_enriched_genes)

    DefaultAssay(seurat_object) <- "RNA"
    dendrogram_avg_expression <- AverageExpression(seurat_object,
                                                   features = final_enriched_genes,
                                                   assays = "RNA")
    dist <- dist(t(dendrogram_avg_expression$RNA), method = "canberra")
    hc <- hclust(dist, method = "complete")
    final_dendrogram <- hc

    # let's use the cleaned up meta data to ensure that celltype classifications
    # match with cluster labels for node label replacement
    nicknames_tibble <- my_final_metadata %>%
      select(final_cluster_labels, celltype) %>%
      dplyr::group_by(final_cluster_labels, celltype) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::mutate(cluster_numeric = as.numeric(final_cluster_labels)) %>%
      dplyr::arrange(cluster_numeric) %>%
      dplyr::mutate(label = paste0(final_cluster_labels, ": ", celltype))

    final_dendrogram$labels <- nicknames_tibble$label
  } else {
    final_dendrogram <- replacement_dendrogram
  }


  #############################################
  #####      save default reduction       #####
  #############################################

  clustering_method <- ifelse(is.null(seurat_object@reductions$umap) == FALSE, "umap",
                       ifelse(is.null(seurat_object@reductions$tsne) == FALSE, "tsne",
                       ifelse(is.null(seurat_object@reductions$pca) == FALSE, "pca", "none")))


  misc <- list(clustering_method = clustering_method,
               final_colors = final_colors,
               final_library_colors = final_library_colors,
               dendrogram = final_dendrogram
               )


  #############################################
  #####   reduce size of the data slot    #####
  #############################################
  # creating a small empty matrix to replace large expression slots
  my_matrix <- matrix(0, 1, 1)
  # the shiny app does not need the scale.data slot
  my_assays$RNA@scale.data <- my_matrix

  # if the object has an integrated slot
  if (is.null(my_assays$integrated) == FALSE) {
    my_assays$integrated@scale.data <- my_matrix
  }

  # the seurat_object_big is used for reclustering, as the integrated
  # data slot is necessary for re-normalization
  seurat_obj_big <- new("Seurat",
                        reductions = my_reductions,
                        meta.data = my_final_metadata,
                        version = my_version,
                        assays = my_assays,
                        active.assay = my_active_assay,
                        active.ident = my_active_ident)

  # if the object has an integrated slot
  if (is.null(my_assays$integrated) == FALSE) {
    my_assays$integrated <- NULL
  }
  # the shiny app (except for the reclustering) does not
  # need the counts slot
  my_assays$RNA@counts <- my_matrix


  #############################################
  #####      Create the Shiny Objects     #####
  #############################################

  seurat_obj <- new("Seurat",
                    reductions = my_reductions,
                    meta.data = my_final_metadata,
                    version = my_version,
                    assays = my_assays,
                    active.assay = my_active_assay,
                    active.ident = my_active_ident,
                    misc = misc)



  #############################################
  #####  save created objects for shiny   #####
  #############################################

  save(seurat_obj, file = paste0(export_data_path, "seurat_obj.RData"))
  save(seurat_obj_big, file = paste0(export_data_path, "seurat_obj_big.RData"))
}
