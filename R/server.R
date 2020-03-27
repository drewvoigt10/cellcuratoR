#' Shiny app server
#'
#' @param input provided by shiny
#' @param output provided by shiny
#'
#' @export
#' @importFrom ggplot2 ggplot aes element_blank theme geom_point geom_text position_nudge
#' @importFrom ggplot2 theme_classic ggtitle xlab ylab geom_tile scale_x_continuous
#' @importFrom ggplot2 scale_fill_gradient scale_y_reverse coord_flip ggsave
#' @importFrom dplyr full_join mutate filter select arrange desc
#' @importFrom DT renderDataTable
#' @importFrom magrittr "%>%"
#' @importFrom grDevices hcl
#' @importFrom stats hclust rnorm
#' @importFrom Seurat Idents DimPlot FeaturePlot GetAssayData NormalizeData FindVariableFeatures WhichCells RunUMAP RunTSNE VariableFeatures ScaleData RunPCA FindNeighbors Tool FindClusters FindMarkers DefaultAssay FetchData HoverLocator AverageExpression
#' @importFrom ggdendro dendro_data segment label
#' @importFrom plotly ggplotly renderPlotly subplot plot_ly layout event_data
#' @importFrom patchwork plot_layout
#' @importFrom shiny Progress validate reactive renderText renderUI req validate conditionalPanel selectInput renderImage actionButton selectizeInput
#' @importFrom shiny updateSelectizeInput downloadButton eventReactive checkboxGroupInput radioButtons sliderInput withProgress shinyServer icon plotOutput
#' @importFrom shiny strong
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom rlang %||% UQ
#' @importFrom shinyFiles shinyDirChoose parseDirPath getVolumes
#' @importFrom methods is new

shinyAppServer <- shinyServer(function(session, input, output) {
  ### Chunk 1: select your dataset and load in associated data objects/RDS files
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  shinyDirChoose(input, "directory", roots = volumes, session = session, restrictions = system.file(package = "base"))

  myProjects <- reactive({
    # Simply list all available projects created with the export seurat function and stored in the data2 directory
    if (is.integer(input$directory)) {
      cat("No directory has been selected (shinyDirChoose)")
    } else {
      list.files(parseDirPath(volumes, input$directory))
    }
  })

  output$directoryWarning <- renderText({
    req(input$directory)
    req(myProjects())
    if(("seurat_obj.RData") %in% myProjects()){
      warn_message <- "WARNING: directory not selected properly. Consider selecting one directory level up."
    } else {
      warn_message <- ""
    }

    warn_message

  })

  output$dataSelect <- renderUI({
    # dataSelect contains the UI allowing the user to select which dataset should be loaded
    req(input$directory)
    conditionalPanel(
      condition = "input.menutab == 'dashboard'",
      selectInput(inputId = "dataset",
                  label = h3("Which dataset should be loaded?"),
                  choices = c("",myProjects()),
                  selected = NULL)
    )


  })

  loaded_data <- reactive({
    # Load the data from the selected dataset (UI = dataSelect, choices = myProjects())
    req(input$dataset)
    req(input$directory)
    load(file.path(as.character(parseDirPath(volumes, input$directory)),  input$dataset, "seurat_obj.RData"))
    seurat_obj
  })

  final_colors <- reactive({
    # Load in the final_colors.rds file created by the export seurat object function
    req(input$dataset)
    final_colors <- loaded_data()@misc$final_colors
    #readRDS(paste0("../data2/", input$dataset, "/", "final_colors.rds"))
    final_colors

  })

  final_library_colors <- reactive({
    # Load in the final_library_colors.rds file created by the export seurat object function
    req(input$dataset)
    final_library_colors <- loaded_data()@misc$final_library_colors
    #readRDS(paste0("../data2/", input$dataset, "/", "final_library_colors.rds"))
    final_library_colors
  })

  hc_final <- reactive({
    # Load in hc_final, the Hierarchical Clustering dendrogram created by the export seurat object function
    req(input$dataset)
    #load(paste0("../data2/", input$dataset, "/", "final_dendrogram.RData"))
    final.dendrogram <- loaded_data()@misc$dendrogram
    final.dendrogram
  })

  all_genes <- reactive({
    # Load in the csv file containing all gene names used in the analysis. Used for quick accessibility
    # of genes in the heatmap and violin plot interactive helpers.
    req(input$dataset)
    genes <- data.frame(rownames(loaded_data()@assays$RNA@data))
    colnames(genes) <- "genes"
    genes
  })

  n_clusters <- reactive({
    # Computes the number of unique clusters in the final_cluster_labels slot of the loaded seurat object (loaded_data())
    n_clusters <-
      loaded_data()$final_cluster_labels %>%
      unique() %>%
      length()
  })

  ### Chunk 2: visualize clustered cells (tabset1)
  output$umapPlot <- renderPlotly({
    # In tabset 1, the umap dimensionality reduction visualization
    validate(
      need(input$dataset != "", "Please select a data set from the left hand menu")
    )

    loaded_plot_data <- loaded_data()
    # Cells are colored according to the selection in the UI tSNE_plot_color
    idents_use <- if(input$tSNE_plot_color == "clusters") {"final_cluster_labels"} else {"libraryID"}
    colors_use <- if(input$tSNE_plot_color == "clusters") {final_colors()} else {final_library_colors()}
    Seurat::Idents(object = loaded_plot_data) <- idents_use
    tsne_plot <- DimPlot(object = loaded_plot_data, cols = colors_use)
    main <- HoverLocator(plot = tsne_plot,
                         information = FetchData(object = loaded_plot_data,
                                                 vars = c("final_cluster_labels", "libraryID", "celltype")))

    # Next we create a legend for the right side of the plot that depicts the color scheme
    legend_label <- levels(loaded_plot_data@active.ident)
    legend_x_cord <- rep(1, length(legend_label))
    legend_y_cord <- rev(seq(1:length(legend_label)))
    manual_legend_data <- data.frame(legend_x_cord, legend_y_cord, legend_label)

    # The legend will have different widths if displaying simply cluster numbers (when input$tSNE_plot_color == "clusters")
    # versus if displaying longer character vectors containing donorID (input$tSNE_plot_color == "libraryID")
    reactive_size <- ifelse(input$tSNE_plot_color == "clusters", 5, 3)
    reactive_width1 <- ifelse(input$tSNE_plot_color == "clusters", 0.9, 0.7)
    reactive_width2 <- ifelse(input$tSNE_plot_color == "clusters", 0.1, 0.3)

    tsne_leg <- ggplot(data = manual_legend_data,
                       mapping = aes(x     = as.factor(legend_x_cord),
                                     y     = legend_y_cord,
                                     label = as.character(legend_label))) +
      geom_point(color = colors_use, size = 5)  +
      geom_text(position = position_nudge(x = 0), size = reactive_size) +
      theme_classic() +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()
      )

    leg <- ggplotly(tsne_leg +
                      theme_classic() +
                      theme(legend.position = "none",
                            axis.line=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks=element_blank()) +
                      xlab("") +
                      ylab(""))
    sp <- subplot(main, leg, widths = c(reactive_width1, reactive_width2), titleY = TRUE, titleX = TRUE)
    sp
  })

  output$umapHelper <- renderUI({
    # The helper UI in which the user can delect how the plot should be colored: according to cell cluster or library
    conditionalPanel(
      condition = "input.tabset1 == 'Dimensionality Reduction'",
      selectInput(inputId = "tSNE_plot_color",
                  label = h4("Dimensionality reduction plot color"),
                  choices = c("clusters", "libraries"),
                  selectize = TRUE,
                  selected = NULL)
    )
  })

  output$featurePlot <- renderPlotly({
    # Returns a feature plot - a heatmap of the dimensionality reduction overlayed with expression of th egene of interest
    req(input$dataset)
    req(input$plot_gene_heatmap)
    loaded_plot_data <- loaded_data()
    Seurat::Idents(object = loaded_plot_data) <- "final_cluster_labels"
    feature_plot = FeaturePlot(object = loaded_plot_data,
                               features = input$plot_gene_heatmap,
                               cols = c("grey", "blue")) +
      theme(legend.position = "none") +
      ggtitle("")

    # Below, we are creating a custom color scale to diplay on the right side of the feature plot
    # this scale depicts the transcripts per 10K (as all data has been normalized with a scale factor of 10,000)
    # This helps in correlating expression in the feature plot to the supplemental spreadsheets created for each publication.
    data <- GetAssayData(object = loaded_plot_data, slot = "data") # Acquire the RNA assay data
    max_expression <- max(expm1(data[input$plot_gene_heatmap, ]))  # Find the cell with the maximum expression of the gene of interest (plot_gene_heatmap)

    if(max_expression > 0) {
      expression_sequenced <- data.frame(TP.10k = 1:round(max_expression))
    } else {
      # If there is no expression of a gene, we still need to create a dataframe that can be converted into a
      # legend so our plots look consistent, regardless of gene expression.
      expression_sequenced <- data.frame(TP.10k = 1:2)
    }

    heatmap_legend <- ggplot(expression_sequenced) +
      geom_tile(aes(x = 1, y = TP.10k, fill = TP.10k)) +
      scale_x_continuous(limits=c(0,2),breaks=1, labels = "TP.10K") +
      theme_classic() +
      theme(legend.position = "none") +
      xlab("") +
      ylab("")

    if(max_expression > 0) {
      heatmap_legend <- heatmap_legend + scale_fill_gradient(low = 'lightgrey', high = 'blue')
    } else {
      heatmap_legend <- heatmap_legend + scale_fill_gradient(low = 'lightgrey', high = 'lightgrey')
    }

    feature_plot <- HoverLocator(plot = feature_plot,
                                 information = FetchData(object = loaded_plot_data,
                                                         vars = c("final_cluster_labels", "libraryID", "celltype")))

    subplot(feature_plot,
            ggplotly(heatmap_legend),
            widths = c(0.9, 0.1),
            titleY = TRUE,
            titleX = TRUE)

  })

  output$featurePlotHelper <- renderUI({
    # The helper UI for the heatmap in which the user can select which gene should be heat-map-ified
    conditionalPanel(
      condition = "input.tabset1 == 'Heatmap'",
      selectizeInput("plot_gene_heatmap",
                     label = h3("gene for heatmap"),
                     choices = NULL,
                     multiple = FALSE,
                     options= list(maxOptions = 50)
      )
    )
  })

  observeEvent(input$dataset, {
    choicesVec <- all_genes()$genes
    updateSelectizeInput(session, "plot_gene_heatmap",
                         choices = choicesVec,
                         selected = NULL,
                         server=TRUE)
  })

  ### Chunk 3: Other visualizations (violin, recluster, DE) (tabset 2)
  my_violin_object <- reactive({
    # the violin plot visualization. The identity class of the seurat object is automatically switched to final cluster labels.
    # the violin figure function takes the seurat object (loaded_data()), gene_to_investigate (input$violin_genes),
    #dendrogram_input (hc_final()), and colors (final_colors()) as inputs
    req(input$violin_genes)
    loaded_plot_data <- loaded_data()
    Seurat::Idents(object = loaded_plot_data) <- "final_cluster_labels"
    genes_to_investigate <- input$violin_genes

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

    p_violins <-
      construct_violin_plot(
        my_object = loaded_plot_data,
        genes_to_investigate = genes_to_investigate,
        dendrogram_input = hc_final(),
        colors = final_colors(),
        use_noise = TRUE
      ) +
      theme_violins()

    dendrogram_data <- ggdendro::dendro_data(hc_final(), type = "rectangle")

    p_dendrogram <-
      ggdendro::segment(dendrogram_data) %>%
      ggplot(aes(x = x, y = y)) +
      geom_segment(
        aes(
          xend = xend,
          yend = yend
        )
      ) +
      geom_text(
        data = ggdendro::label(dendrogram_data),
        aes(label = label),
        hjust = 1,
        vjust = -0.7,
        size = 6
      ) +
      scale_y_reverse() +
      scale_x_continuous(expand = c(0, 0.6)) +
      coord_flip() +
      theme_violins()

    p_dendrogram +
      p_violins +
      patchwork::plot_layout(ncol = 2, widths = c(3, 7))
  })

  output$violin_plot <- renderPlot({
    # the violin plot visualization. The identity class of the seurat object is automatically switched to final cluster labels.
    # the violin figure function takes the seurat object (loaded_data()), gene_to_investigate (input$violin_genes),
    #dendrogram_input (hc_final()), and colors (final_colors()) as inputs
    my_violin_object()
  })

  output$violin_download <- downloadHandler(
    filename = function(){paste("violin_plot",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=my_violin_object(), height = 9, width = 2*length(input$violin_genes) + 1, units = "in")
    }
  )

  output$violinHelper <- renderUI({
    # The helper UI for the violin plot. Multiple genes can be input if the violin plot tab is selected within tabset 2.
    conditionalPanel(
      condition = "input.tabset2 == 'Violin Plot'",
      selectizeInput('violin_genes',
                     label = h3("Violin genes"),
                     choices = NULL,
                     multiple=TRUE,
                     options= list(maxOptions = 50)
      )

    )

  })

  observeEvent(input$dataset, {
    choicesVec <- all_genes()$genes
    updateSelectizeInput(session, "violin_genes",
                         choices = choicesVec,
                         selected = NULL,
                         server=TRUE)
  })

  output$violinDownload <- renderUI({
    # Download handler for the violin plot
    conditionalPanel(
      condition = "input.tabset2 == 'Violin Plot'",
      downloadButton(outputId = "violin_download",
                     label = "Download the plot")
    )
  })

  ######################
  ### RECLUSTER CODE ###
  ######################
  loaded_data_big <- eventReactive(input$start_recluster, {
    # a larger seurat object is required for reclustering, due to the requirement that reclustering accesses the integrated data slot.
    # In order to keep the app fast, this larger seurat object (which is automatically created by the export seurat function) is
    # only loaded in if the user requests to recluster the data.
    req(input$dataset)
    load(file.path(as.character(parseDirPath(volumes, input$directory)), input$dataset, "seurat_obj_big.RData"))
    seurat_obj_big
  })

  recluster_object <- eventReactive(input$start_recluster, {
    # using the seurat_recluster function, takes the larger seurat object (loaded_data_big()) and
    # Pulls out cells belonging to subset_clusters and performs renormalization followed by
    # UMAP, tSNE, and PCA dimensionality reduction.
    req(loaded_data_big())

    # Create a Progress object
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Reclustering", value = 0)

    seurat_object <- loaded_data_big()
    Seurat::Idents(seurat_object) <- "final_cluster_labels"

    new_object <- subset(seurat_object, idents = as.character(input$select_recluster))

    ## Re-normalize (UMAP requires a large seurat object with the integrated data slot)
    new_object <- NormalizeData(new_object, normalization.method = "LogNormalize", scale.factor = 10000)
    progress$inc(1/6, detail = paste("Finding Variable Features"))

    new_object <- FindVariableFeatures(new_object, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(new_object)
    progress$inc(1/6, detail = paste("Scaling Data"))
    new_object <- ScaleData(new_object, features = all.genes)

    progress$inc(1/6, detail = paste("Running PCA"))
    new_object <- RunPCA(new_object, features = VariableFeatures(object = new_object))

    progress$inc(1/6, detail = paste("Finding Clusters and Neighbors"))
    new_object <- FindNeighbors(new_object, dims = 1:20)
    new_object <- FindClusters(new_object, resolution = 0.5)

    progress$inc(1/6, detail = paste("Running UMAP"))
    new_object <- RunUMAP(new_object, dims = 1:20)

    progress$inc(1/6, detail = paste("Running tSNE"))
    new_object <- RunTSNE(new_object, dims = 1:20)


    Seurat::Idents(new_object) <- "final_cluster_labels"

    return(new_object)
  })

  output$recluster_plot <- renderPlotly({
    # The plot to visualize the reclustered data. Cells in the reclustered plot can be visualized by cluster, libraryID, or gene,
    # which is contained in the input$reductionPlotColor ui-side input.
    # The dimensionality reduction used can be modified by changing the input$recluster_reduction switch.
    validate(
      need(input$recluster_reduction != "", "Please select a dimensionality reduction technique for reclustering (this can be changed later without requiring re-computation)")
    )
    req(recluster_object())
    recluster_object <- recluster_object()

    if(input$reductionPlotColor == 1){ ### color cells by CLUSTER
      Seurat::Idents(recluster_object) <- "final_cluster_labels"
      label_index <- names(summary(loaded_data()@meta.data$final_cluster_labels))
      selected_groups_index <- as.character(input$select_recluster)
      p <- DimPlot(recluster_object,
                   reduction = input$recluster_reduction,
                   cols = final_colors()[which(label_index %in% selected_groups_index)])
    }

    if(input$reductionPlotColor == 2){ ### color cells by LIBRARYID
      Seurat::Idents(recluster_object) <- "libraryID"
      p <- DimPlot(recluster_object,
                   reduction = input$recluster_reduction,
                   cols = final_library_colors())
    }

    if(input$reductionPlotColor == 3){ ### color cells by GENE
      req(input$reduction_plot_gene_heatmap)
      Seurat::Idents(recluster_object) <- "final_cluster_labels"
      feature_plot <- FeaturePlot(object = recluster_object,
                                  features = input$reduction_plot_gene_heatmap,
                                  reduction = input$recluster_reduction,
                                  cols = c("grey", "blue")) +
        theme(legend.position = "none") +
        ggtitle("")

      # Below, we are creating a custom color scale to diplay on the right side of the feature plot
      # this scale depicts the transcripts per 10K (as all data has been normalized with a scale factor of 10,000)
      # This helps in correlating expression in the feature plot to the supplemental spreadsheets created for each publication.
      data <- GetAssayData(object = recluster_object, slot = "data") # Acquire the RNA assay data
      max_expression <- max(expm1(data[input$reduction_plot_gene_heatmap, ]))  # Find the cell with the maximum expression of the gene of interest (plot_gene_heatmap)

      if(max_expression > 0) {
        expression_sequenced <- data.frame(TP.10k = 1:round(max_expression))
      } else {
        # If there is no expression of a gene, we still need to create a dataframe that can be converted into a
        # legend so our plots look consistent, regardless of gene expression.
        expression_sequenced <- data.frame(TP.10k = 1:2)
      }

      heatmap_legend <- ggplot(expression_sequenced) +
        geom_tile(aes(x = 1, y = TP.10k, fill = TP.10k)) +
        scale_x_continuous(limits=c(0,2),breaks=1, labels = "TP.10K") +
        theme_classic() +
        theme(legend.position = "none") +
        xlab("") +
        ylab("")

      if(max_expression > 0) {
        heatmap_legend <- heatmap_legend + scale_fill_gradient(low = 'lightgrey', high = 'blue')
      } else {
        heatmap_legend <- heatmap_legend + scale_fill_gradient(low = 'lightgrey', high = 'lightgrey')
      }

      feature_plot <- HoverLocator(plot = feature_plot,
                                   information = FetchData(object = recluster_object,
                                                           vars = c("final_cluster_labels", "libraryID", "celltype")))

      p <- subplot(feature_plot,
                   ggplotly(heatmap_legend),
                   widths = c(0.9, 0.1),
                   titleY = TRUE,
                   titleX = TRUE)
    }
    p

  })

  list_of_groups <- reactive({
    # Reactive expression that lists all of the final_cluster_labels within the loaded_data() seurat object
    loaded_data()@meta.data$final_cluster_labels %>%
      levels()
  })

  output$reclusterHelper <- renderUI({
    # The UI helper for reclustering. This helper displays list_of_groups(), which lists all of the
    # final_cluster_labels in the loaded seurat object.
    conditionalPanel(
      condition = "input.tabset2 == 'Recluster'",
      checkboxGroupInput(inputId = "select_recluster",
                         label = h3("select groups to recluster"),
                         choices = list_of_groups(),
                         selected = NULL)
    )
  })

  output$reclusterHelper2 <- renderUI({
    # The second UI helper for reclustering. This helper displays a reactive start button that the user can
    # select once the groups for reclustering (input$select_recluster) and the dimensionality reduction strategy
    # (input$recluster_reduction) have been selected.

    conditionalPanel(
      condition = "input.tabset2 == 'Recluster'",
      actionButton(inputId = "start_recluster",
                   label = "Start reclustering (slow)")
    )
  })

  output$reclusterHelper3 <- renderUI({
    # The third UI helper for reclustering. This helper displays the choices for dimensionality reduction.
    conditionalPanel(
      condition = "input.tabset2 == 'Recluster'",
      radioButtons("recluster_reduction",
                   label = "Dimensionality Reduction for reclustering",
                   choices = c("pca", "tsne", "umap"),
                   selected = character(0)))
  })

  output$reclusterHelper4 <- renderUI({
    # The fourth UI helper for reclustering. This helper displays the options for how the reclustered cells
    # should be colored (by original cluster, by libraryID, or a heatmap of the gene).
    conditionalPanel(
      condition = "input.tabset2 == 'Recluster'",
      selectInput(inputId = "reductionPlotColor",
                  label = h3("How should the cells be colored?"),
                  choices = list("by cluster" = 1, "by libraryID" = 2, "by gene" = 3),
                  selected = 1)
    )
  })

  output$reclusterHelper5 <- renderUI({
    # The fifth UI helper for reclustering. If the user selects that reclustered cells should be colored by
    # the expression of a particular gene (input$reductionPlotColor == 3), then another UI pops up
    # where the user can input the gene of interest.

    conditionalPanel(
      condition = "input.tabset2 == 'Recluster' && input.reductionPlotColor == 3",
      selectizeInput("reduction_plot_gene_heatmap",
                     label = h3("gene for heatmap"),
                     choices = NULL,
                     multiple = FALSE,
                     selected = NULL,
                     options= list(maxOptions = 50)
      )
    )
  })

  observeEvent(input$dataset, {
    choicesVec <- all_genes()$genes
    updateSelectizeInput(session, "reduction_plot_gene_heatmap",
                         choices = choicesVec,
                         selected = NULL,
                         server=TRUE)
  })

  ######################
  ### DE CODE ###
  ######################
  plot_data <- reactive({
    # plot_data() stores the UMAP or tSNE coordinates of genes for interactive selection with plotly functions.
    # For example, the user can select the laso tool and circle populations of interest to perform differential expression.
    # In order to streamline the selection and subsequent differential expression of selected cells, relevant information
    # including the barcode of selected cells, the final_cluster_labels, and other meta data is appended to plot data
    # which is then used to create the relevant differential-expression plots (dePlot1 and dePlot2, below).
    if(input$tabset1 == "Differential Expression (Group 1)" | input$tabset2 == "Differential Expression (Group 2)"){
      loaded_data <- loaded_data()
      clustering_method <- loaded_data@misc$clustering_method ## default clustering method used,
      ## either "umap" or "tsne" - automatically generated by the export seurat object function

      if(clustering_method == "umap"){
        barcode = rownames(loaded_data@reductions$umap@cell.embeddings)
        embeddings = loaded_data@reductions$umap@cell.embeddings
        plot_data = data.frame(barcode, embeddings)
      } else { ##clustering_method == "tsne"
        barcode = rownames(loaded_data@reductions$tsne@cell.embeddings)
        embeddings = loaded_data@reductions$tsne@cell.embeddings
        plot_data = data.frame(barcode, embeddings)
      }

      ## append on meta data
      barcode <- rownames(loaded_data@meta.data)
      relevant_metadata <- loaded_data@meta.data[,c("final_cluster_labels", "libraryID", "celltype")]
      plot_metadata <- data.frame(barcode, relevant_metadata)
      colnames(plot_metadata) <- c("barcode", "final_cluster_labels", "libraryID", "celltype")

      plot_data <- plot_data %>%
        dplyr::full_join(plot_metadata) %>%
        tibble::as_tibble()

      plot_data
    }
  })

  output$dePlot1 <- renderPlotly({
    # Dimensionality reduction plot displayed in tabset 1. Used for selection of differential expression groups.

    validate(
      need(input$tabset2 == "Differential Expression (Group 2)", "Please select the differential expression tab in the box to the right"),
      need(input$DE_class != "", "Please select differential expression criteria below, select clusters, and hit the \"start differential expression analysis\" button.")
    )

    req(input$DE_class)
    clustering_method <- loaded_data()@misc$clustering_method

    x_axis <- list(
      title = paste0(clustering_method, "_1"),
      zeroline = FALSE
    )

    y_axis <- list(
      title = paste0(clustering_method, "_2"),
      zeroline = FALSE
    )

    plot_vars <- colnames(plot_data())

    p <- plot_ly(plot_data(),
                 x = as.data.frame(plot_data())[,2], # Dimensionality reduction coordinate 1
                 y = as.data.frame(plot_data())[,3], # Dimensionality reduction coordinate 2
                 type = "scatter",
                 mode = "markers",
                 color = ~final_cluster_labels,
                 colors = final_colors(),
                 hoverinfo = 'text',
                 hovertext = paste("cluster :", plot_data()$final_cluster_labels,
                                   "<br> class :", plot_data()$celltype),
                 source = "plot_group_1",
                 key= ~barcode) %>%
      layout(xaxis=x_axis, yaxis=y_axis, showlegend = FALSE)

    if(input$DE_class == 2) {
      p <- p %>% layout(dragmode = "lasso")
    } ## select the lasso tool if the user selects "draw lasso"
    p

  })

  output$dePlot2 <- renderPlotly({
    ## Dimensionality recution plot displayed in tabset 2. Used for selection of differential expression groups.
    req(plot_data())
    req(input$dataset)
    req(input$DE_class)

    if(input$DE_class == 2 & input$DE_identity == 1){
      clustering_method <-  loaded_data()@misc$clustering_method

      x_axis <- list(
        title = paste0(clustering_method, "_1"),
        zeroline = FALSE
      )

      y_axis <- list(
        title = paste0(clustering_method, "_2"),
        zeroline = FALSE
      )

      p <- plot_ly(plot_data(),
                   x = as.data.frame(plot_data())[,2],
                   y = as.data.frame(plot_data())[,3],
                   type = "scatter",
                   mode = "markers",
                   color = ~final_cluster_labels,
                   colors = final_colors(),
                   hoverinfo = 'text',
                   hovertext = paste("cluster :", plot_data()$final_cluster_labels,
                                     "<br> class :", plot_data()$celltype),
                   source = "plot_group_2", key= ~barcode) %>%
        layout(xaxis=x_axis, yaxis=y_axis, showlegend = FALSE)

      if(input$DE_class == 2) {
        p <- p %>% layout(dragmode = "lasso")
      } ## select the lasso tool if the user selects "draw lasso"

      p
    }
  })

  lasso_selection_1 <- reactive({
    # Contains the barcodes of cells selected in lasso 1 in dePlot1
    lasso_1_data <- event_data(event = "plotly_selected", source = "plot_group_1")
    plot_1_filter <- plot_data()[plot_data()$barcode %in% lasso_1_data$key,]
    as.character(plot_1_filter$barcode)
  })

  lasso_selection_2 <- reactive({
    # Contains the barcodes of cells selected in lasso 2 in dePlot2
    lasso_2_data <- event_data(event = "plotly_selected", source = "plot_group_2")
    plot_2_filter <- plot_data()[plot_data()$barcode %in% lasso_2_data$key,]
    as.character(plot_2_filter$barcode)
  })

  output$deHelper <- renderUI({
    # The first helper function displayed under tabset 1. Contains the logFC threshold for the
    # FindMarkers function (which performs the DGE). Higher logFC threshold results in filtering of smaller
    # gene expression differences and runs faster (but is obviously less sensitive).
    conditionalPanel(
      condition = "input.tabset1 == 'Differential Expression (Group 1)'",
      sliderInput(inputId = "logFC_threshold",
                  label = "logFC threshold",
                  min = 0,
                  max = 1,
                  step = 0.05,
                  value = 0.50)
    )
  })

  output$deHelper2 <- renderUI({
    # The second helper function displayed under tabset 1. Contains the min.pct threshold for the
    # FindMarkers function (which performs the DGE). Higher min.pct threshold results in filtering of
    # genes that are expressed in a low percentage of cells and runs faster (but is obviously less sensitive).
    conditionalPanel(
      condition = "input.tabset1 == 'Differential Expression (Group 1)'",
      sliderInput(inputId = "min_pct_threshold",
                  label = "minimum percent threshold",
                  min = 0,
                  max = 1,
                  step = 0.05,
                  value = 0.50)
    )
  })

  output$deHelper3 <- renderUI({
    # The third helper function displayed under tabset 1. Allows the user to select the comparison
    # groups for differential expression. The user can perform DGE on clusters of interest, for example
    # clusters 1,4,5 vs clusters 2,7, 9. Or, the user can use the lasso tool (2) to manually circle cells
    # for comparisons.

    conditionalPanel(
      condition = "input.tabset1 == 'Differential Expression (Group 1)'",
      radioButtons(inputId = "DE_class",
                   label = "Select cells for differential expression:",
                   choices = c("predefined clusters" = 1,
                               "manually selected populations (draw laso)" = 2),
                   selected = character(0))
    )
  })

  output$deHelper5 <- renderUI({
    # The fifth helper function, which is displayed under tabset 1. Allows the user to perform differential
    # expression between different classes of cells in the meta.data. For example, differential expression
    # can be performed between cells originating from the fovea vs periphery if a binary "region" variable
    # is exported with the meta.data using the export_seurat_object() function.

    additional_idents <- colnames(loaded_data()@meta.data)[-c(1:3)]
    conditionalPanel(
      condition = "input.tabset1 == 'Differential Expression (Group 1)'",
      radioButtons(inputId = "DE_identity",
                   label = "Perform differential expression between:",
                   choices = c("all selected cells (default)" = 1,
                               additional_idents),
                   selected = 1)
    )
  })

  output$DE_group_1 <- renderUI({
    # The first helper function displayed under tabset 2. If input$DE_class == 1 (the user wants to perform
    # DGE on predefined clusters of cells), then a helper appears under tabset 2 listing the comparison groups
    # in a checkbox format. Groups selected here form the comparison group "1" in the FindMarkers function.
    # Of note, list_of_groups() lists all of the final_cluster_labels in the loaded seurat object.

    conditionalPanel(
      condition = "input.tabset2 == 'Differential Expression (Group 2)' && input.DE_class == 1",
      checkboxGroupInput(inputId = "DE_class_1",
                         label = h3("Group 1"),
                         choices = list_of_groups(),
                         selected = NULL)
    )
  })

  output$DE_group_2 <- renderUI({
    # The second helper function displayed under tabset 2. If input$DE_class == 1 (the user wants to perform
    # DGE on predefined clusters of cells), then a helper appears under tabset 2 listing the comparison groups
    # in a checkbox format. Groups selected here form the comparison group "2" in the FindMarkers function.
    # Of note, list_of_groups() lists all of the final_cluster_labels in the loaded seurat object.
    # If the user selects to perform differential expression on an additional identity (ie region: input$DE_identity != 1),
    # then this panel will NOT be displayed.
    conditionalPanel(
      condition = "input.tabset2 == 'Differential Expression (Group 2)' && input.DE_class == 1 && input.DE_identity == 1",
      checkboxGroupInput(inputId = "DE_class_2",
                         label = h3("Group 2"),
                         choices = list_of_groups(),
                         selected = NULL)
    )
  })

  output$deHelper4 <- renderUI({
    ## The fourth helper function - contains an action button to start differential expression
    conditionalPanel(
      condition = "input.tabset2 == 'Differential Expression (Group 2)'",
      actionButton("start_DE", "Start Differential Expression Analysis")
    )
  })

  DE_data <- eventReactive(input$start_DE, {
    # DE_data contains an additional column in the meta.data, de_group, that communicates what cells
    # should be in each investigated group for differential expression.

    loaded_data <- loaded_data()
    if(input$DE_class == 1) { # If the user wants to perform DGE on predefined clusters of cells...

      if(input$DE_identity == 1){ # If the user wants to perform DGE between cell selections...
        Idents(loaded_data) <- "final_cluster_labels"
        de_group <- ifelse(loaded_data@meta.data$final_cluster_labels %in% input$DE_class_1, "Group_1",
                           ifelse(loaded_data@meta.data$final_cluster_labels %in% input$DE_class_2, "Group_2", "none"))
      } else { # If the user wants to perform DGE between a binary biological variable...
        # Although it would be easier to subset loaded_data to all cells in input$DE_class_1,
        # subset requires an integrated data slot (and would necessitate loading the large seurat
        # object, which I don't want to do here). I therefore implement this work-around:
        Idents(loaded_data) <- as.character(input$DE_identity)
        binary_idents = levels(loaded_data@active.ident)

        de_group <- loaded_data@meta.data %>%
          mutate(de_group = ifelse((final_cluster_labels %in% input$DE_class_1) & (UQ(as.symbol(as.character(input$DE_identity))) == binary_idents[1]), "Group_1",
                                   ifelse((final_cluster_labels %in% input$DE_class_1) & (UQ(as.symbol(as.character(input$DE_identity))) == binary_idents[2]), "Group_2", "none"))
          ) %>%
          dplyr::select(de_group)
      }

      loaded_data@meta.data <- data.frame(loaded_data@meta.data, de_group)

    } else if(input$DE_class == 2) { # If the user wants to perform DGE on lasso selections...
      if(input$DE_identity == 1){ # If the user wants to perform DGE between cell selections...
        req(lasso_selection_1())
        req(lasso_selection_2())
        Idents(loaded_data) <- "final_cluster_labels"
        barcodes <- rownames(loaded_data@meta.data)
        de_group <- ifelse(barcodes %in% lasso_selection_1(), "Group_1",
                           ifelse(barcodes %in% lasso_selection_2(), "Group_2", "not_selected"))

      } else { # If the user wants to perform DGE between a binary biological variable...
        req(lasso_selection_1())
        barcodes <- rownames(loaded_data@meta.data)
        Idents(loaded_data) <- as.character(input$DE_identity)
        binary_idents = levels(loaded_data@active.ident)

        de_group <- loaded_data@meta.data %>%
          rownames_to_column("barcodes") %>%
          as_tibble() %>%
          mutate(de_group = ifelse((barcodes %in% lasso_selection_1()) & (UQ(as.symbol(as.character(input$DE_identity))) == binary_idents[1]), "Group_1",
                                   ifelse((barcodes %in% lasso_selection_1()) & (UQ(as.symbol(as.character(input$DE_identity))) == binary_idents[2]), "Group_2", "not_selected"))) %>%
          dplyr::select(de_group)
      }

      loaded_data@meta.data <- data.frame(loaded_data@meta.data, barcodes, de_group)

    }

    loaded_data@meta.data$de_group <- factor(loaded_data@meta.data$de_group, levels = c("Group_1", "Group_2", "not_selected"))
    Seurat::Idents(loaded_data) <- "de_group"
    loaded_data
  })

  DE_markers <- reactive({
    loaded_data = DE_data()
    Idents(loaded_data) <- "de_group"
    withProgress(message = "performing DGE", value = 0, {
      variables <- find_dge_variables(object = loaded_data,
                                      ident.1 = "Group_1",
                                      ident.2 = "Group_2",
                                      my_logfc_threshold = input$logFC_threshold,
                                      my_minpct_threshold = input$min_pct_threshold)

      chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) # https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r

      if(length(variables) > 20){
        n <- 20
      } else {
        n <- length(variables)
      }

      variable_chunks <- chunk(variables, n)

      markers <- data.frame()
      for(i in 1:n){
        dge.chunk <- FindMarkers(object = loaded_data, ident.1 = "Group_1", ident.2 = "Group_2",
                                 my_logfc_threshold = input$logFC_threshold, my_minpct_threshold = input$min_pct_threshold,
                                 features = variable_chunks[[i]]) %>% tibble::rownames_to_column("gene")
        markers <- rbind(markers, dge.chunk)
        incProgress(1/n, detail = paste(i*(100/n), "% of genes"))
      }

    }
    )

    markers

  })

  output$DE_table <- DT::renderDataTable(
    # Displays differential expression results from DE_markers() in a DT
    DE_markers(), server = FALSE, extensions = c("Buttons"), options = list(scrollX = TRUE, dom = 'Bfrtip', buttons = c('csv', 'excel'))
  )

  output$delta_plotly <- renderPlotly({
    # Displays a delta_plot of differential expression results:
    # x-axis: Delta-percent (pct cells expressing gene in group 1 - pct cells expressing gene in group 2)
    # y-axis: logFC (+ = increased in group 1, - = increased in group 2)
    regional_differences <- DE_markers() %>%
      mutate(delta_percent = pct.1 - pct.2) %>%
      mutate(group = ifelse(avg_logFC > 0, "Group 1", "Group 2"))

    x_axis_delta <- list(
      title = "delta percent",
      tickvals = c(-1,-0.5, 0, 0.5, 1)
    )

    y_axis_delta <- list(
      title = "avg_logFC",
      tickvals = c(-3,-2,-1,0,1,2,3)
    )

    delta_plot <- plot_ly(type = "scatter", data = regional_differences, x = regional_differences$delta_percent, y = regional_differences$avg_logFC,
                          mode = "markers",
                          color = ~group, colors = c("#F8766D", "#00BFC4"),
                          marker = list(size = 10,
                                        line = list(color = "black", width=1)),
                          text = ~paste("gene: ", gene,
                                        '</br> logFC: ', round(avg_logFC, 3),
                                        '</br> delta percent: ', round(delta_percent, 3))) %>% layout(xaxis = x_axis_delta, yaxis = y_axis_delta) %>%
      layout(showlegend = FALSE)
  })

  output$cluster_visualization <- renderPlotly({
    # Displays a dimensionality reduction plot where group 1 is colored in red while group 2 is colored in blue
    req(DE_data())
    loaded_data <- DE_data()

    if(input$DE_identity == 1){ # If the user wants to perform DGE between a selected populations of cells...
      Idents(loaded_data) <- "de_group"
      myplot <- DimPlot(object = loaded_data,
                        cols = c("#F8766D", "#00BFC4", "darkgrey"))
    } else {
      Idents(loaded_data) <- as.character(input$DE_identity)
      binary_idents <- levels(loaded_data@active.ident)
      plot_label <- ifelse(loaded_data$de_group == "Group_1", binary_idents[1],
                           ifelse(loaded_data$de_group == "Group_2", binary_idents[2], 'not_selected'))
      loaded_data@meta.data <- data.frame(loaded_data@meta.data, plot_label)
      loaded_data@meta.data$plot_label <- factor(loaded_data@meta.data$plot_label,
                                                 levels = c(binary_idents[1], binary_idents[2], "not_selected"))
      Idents(loaded_data) <- "plot_label"
      myplot <- DimPlot(object = loaded_data,
                        cols = c("#F8766D", "#00BFC4", "darkgrey"))
    }
    myplot

  })

  #######################
  ### Embedded images ###
  #######################
  output$cellcuratoR <- renderImage({
    return(list(
      src = "../www/cellcuratoR.png",
      filetype = "image/png",
      width = 346.8, #2312
      height = 401.0, #2673
      alt = "error: www directory not found"
    ))
  }, deleteFile = FALSE)



})
