#' Prepare data for violin plots
#'
#' \code{prepare_violin_data_colors} returns a list of expression data for the
#' requested genes as well as an ordered vector of colors that will be used to
#' fill the violin plots. This data can then be plotted by the
#' \code{\link{construct_violin_plot}} function.
#'
#' @param my_object S4 object from which expression data should be extracted.
#'
#' @param genes_to_investigate a vector containing character strings of the
#'   genes that the user would like to visualize with violin plots.
#'
#' @param dendrogram_input an input dendrogram that relates the expression
#'   profiles in each cluster. See \code{\link{export_shiny_object}} for an
#'   example of generating such a dendrogram.
#'
#' @param colors a vector of length n_clusters to color the violin plot. The
#'   order of colors corresponds to the order of final_cluster_labels.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' prepare_violin_data_colors(my_object = seurat_obj, genes_to_investigate = c("RHO", "PDE6A"), dendrogram_input = my_dendrogram, colors = c("red", "green", "orange"))
#'}
#'
#' @export
prepare_violin_data_colors <-
  function(my_object,
           genes_to_investigate,
           dendrogram_input,
           colors) {
    n_genes <- length(genes_to_investigate)

    n_clusters <-
      my_object$final_cluster_labels %>%
      unique() %>%
      length()

    if(class(dendrogram_input) == "hclust"){
      # If the user supplies a dendrogram generated from the export_shiny_object() function,
      # reorder the colors according to the order of leaves in the dendrogram.
      if(grepl(":", dendrogram_input$labels[1], fixed = TRUE)) {
        # If the dendrogram has a colon in the label, then the labels were created by the export seurat function
        # we strip off the colon to save text space when the dendrogram will be plotted.
        ordered_colors <- colors[dendrogram_input$order]

        new_levels <-
          (sub(":.*", "", dendrogram_input$labels[dendrogram_input$order]))

      } else {
        # If labels dont have a colon (:), then they are already ordered by the user.
        ordered_colors <- rev(colors)

        new_levels <-
          Seurat::Idents(object = my_object) %>%
          levels() %>%
          rev()
      }
    } else {
      ordered_colors <- rev(colors)

      new_levels <-
        Seurat::Idents(object = my_object) %>%
        levels() %>%
        rev()
    }

    ordered_colors <-
      purrr::set_names(
        ordered_colors,
        as.character(new_levels)
      )

    Seurat::Idents(object = my_object) <-
      Seurat::Idents(object = my_object) %>%
      factor(levels = new_levels)

    cell_clusters <-
      tibble::tibble(
        cluster_id = Seurat::Idents(my_object),
        cell_id = names(Seurat::Idents(my_object))
      )

    fetch_data <-
      Seurat::FetchData(
        object = my_object,
        vars = genes_to_investigate,
        slot = "data"
      ) %>%
      tibble::rownames_to_column("cell_id")

    # Extracted the noise code from `Seurat:::SingleExIPlot`
    # (cf. lines 18-31). By setting the random number generator seed,
    # Seurat adds a tiny bit of noise to the expression values in the
    # plot. The same noise vector gets added to each gene. I've adapted
    # the strategy below to achieve the same plots as Seurat. The noise
    # changes how the violins get rendered.
    seed.use <- 42

    set.seed(seed = seed.use)

    noise <- rnorm(n = nrow(fetch_data)) / 1e+05

    seurat_data <-
      fetch_data %>%
      tidyr::gather(gene_id, expression_value, -cell_id) %>%
      # Maintain the ordering of the genes as requested in the vector.
      # Otherwise, the genes will be ordered alphabetically.
      dplyr::mutate(gene_id = factor(gene_id, levels = genes_to_investigate)) %>%
      dplyr::inner_join(cell_clusters) %>%
      # Add the noise vector generated above.
      dplyr::mutate(expression_value_noise = expression_value + rep(noise, n_genes))

    list(
      seurat_data = seurat_data,
      ordered_colors = ordered_colors
    )
  }


#' Generate violin plots of prepared data
#'
#' \code{construct_violin_plot} Returns a ggplot object of expression profiles
#' for all clusters according to the input genes of interest. Colors and
#' expression profiles of data used for this visualization are prepared in
#' \code{\link{prepare_violin_data_colors}}.
#'
#' @param my_object S4 object from which expression data should be extracted.
#'   Used in generating data in \code{\link{prepare_violin_data_colors}}.
#'
#' @param genes_to_investigate a vector containing character strings of the
#'   genes that the user would like to visualize with violin plots. Used in
#'   generating data in \code{\link{prepare_violin_data_colors}}.
#'
#' @param dendrogram_input an input dendrogram that relates the expression
#'   profiles in each cluster. See \code{\link{export_shiny_object}} for example
#'   of generating such a dendrogram. Used in generating data in
#'   \code{\link{prepare_violin_data_colors}}.
#'
#' @param colors an vector of length n_clusters to color the violin plot. The
#'   order of colors corresponds to the order of final_cluster_labels. Used in
#'   generating data in \code{\link{prepare_violin_data_colors}}.
#'
#' @param use_noise adds random noise to expression values if TRUE (default).
#'   This results in violins being drawn for only clusters that contain at least
#'   25% of cells expressing the gene of interest above background. If FALSE,
#'   uses unmodified expression profiles. This results in violins being drawn
#'   for all clusters, regardless of the proportion of cells that express the
#'   gene
#'
#' @import ggplot2
#'
#' @param scale the scale for facet_grid
#'
#' @examples
#' \dontrun{
#' construct_violin_plot(my_object = seurat_obj, genes_to_investigate = c("RHO", "PDE6A"),
#' dendrogram_input = my_dendrogram, colors = c("red", "green", "orange"))
#'}
#'
#' @export
construct_violin_plot <- function(my_object,
                                  genes_to_investigate,
                                  dendrogram_input,
                                  colors,
                                  use_noise = TRUE,
                                  scale = "free_x") {

  object_data <-
    prepare_violin_data_colors(
      my_object = my_object,
      genes_to_investigate = genes_to_investigate,
      dendrogram_input = dendrogram_input,
      colors = colors
    )

  p_cluster_expression <-
    object_data$seurat_data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = cluster_id,
        y = if(use_noise) expression_value_noise else expression_value
      )
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = cluster_id),
      scale = "width",
      adjust = 1
    ) +
    ggplot2::facet_grid(~ gene_id, scale = scale) +
    ggplot2::scale_fill_manual(
      values = object_data$ordered_colors
    ) +
    ggplot2::coord_flip()

  p_cluster_expression
}


theme_violins <- function() {
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = rel(1.5)),
    legend.position = "none"
  )
}


