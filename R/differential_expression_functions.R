# Function to get all the descendants on a tree left of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants left of the given node
#
GetLeftDescendants <- function(tree, node) {
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  if (daughters[1] <= (tree$Nnode + 1)) {
    return(daughters[1])
  }
  daughter.use <- GetDescendants(tree, daughters[1])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

# Function to get all the descendants on a tree right of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants right of the given node
#
GetRightDescendants <- function(tree, node) {
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  if (daughters[2] <= (tree$Nnode + 1)) {
    return(daughters[2])
  }
  daughter.use <- GetDescendants(tree = tree, node = daughters[2])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}


#' Find genes that meet selected differential expression thresholds
#'
#' \code{find_dge_variables} returns a vector of genes that satisfy differential
#' expression threshold criteria (log fold-change and min.pct thresholds).
#' This allows for differential expression to be conducted in a loop and a
#' progress bar to display within the shiny app indiciating the progress of the
#' (sometimes lengthy) differential expression. It is modified from the Seurat
#' FindMarkers() function.
#'
#' @param object Seurat S4 object. The active identity of the object should
#' match the ident.1 and ident.2 parameters. You can switch the identity of
#' your Seurat object with
#' \code{Idents(object) <- my_ident}.
#'
#' @param ident.1 Identity class for group-1 for differential expression.
#'
#' @param ident.2 Optional identity class for group-2 for differential
#' expression.If differential expression should be compared to all other
#' cells, leave ident.2 empty.
#'
#' @param my_logfc_threshold The log fold-change threshold that will limit
#' differential expression testing to genes with expression differences
#' greater than the specified log fold-change threshold. Only genes with a
#' log fold-change threshold of greater than \code{my_logfc_threshold} will
#' be returned.
#'
#' @param my_minpct_threshold The cutoff for the proportion of cells that
#' express the gene of interest. Only genes that are expressed in a greater
#' proportion of cells than the \code{my_minpct_threshold} will be returned.
#'
#' @examples
#' \dontrun{
#' find_dge_variables(seurat_obj, ident.1 = 2, ident.2 = c(3,4),
#' my_logfc_threshold = 0.25, my_minpct_threshold = 0.5)
#' find_dge_variables(seurat_obj, ident.1 = 3, my_logfc_threshold = 0.5,
#' my_minpct_threshold = 0.1)
#' }

find_dge_variables <- function(object,
                               ident.1,
                               ident.2,
                               my_logfc_threshold,
                               my_minpct_threshold) {
  group.by = NULL
  assay = NULL
  test.use = "wilcox"
  slot = 'data'
  reduction = NULL
  latent.vars = NULL
  features = NULL

  data.slot <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = 'counts',
    no = slot
  )

  assay <- assay %||% DefaultAssay(object = object)
  data.use <-  GetAssayData(object = object[[assay]], slot = data.slot)

  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(ident.1, ident.2)
    )
  }
  counts <- switch(
    EXPR = data.slot,
    'scale.data' = GetAssayData(object = object[[assay]], slot = "counts"),
    numeric()
  )

  ### FindMarkers.default
  object = data.use
  slot = data.slot
  counts = counts
  cells.1 = ident.1
  cells.2 = ident.2
  features = features
  reduction = reduction
  logfc.threshold = my_logfc_threshold
  test.use = test.use
  min.pct = my_minpct_threshold
  min.diff.pct = -Inf
  verbose = TRUE
  only.pos = FALSE
  max.cells.per.ident = Inf
  random.seed = 1
  latent.vars = NULL
  min.cells.group = 3
  pseudocount.use = 1


  features <- features %||% rownames(x = object)
  methods.noprefiliter <- c("DESeq2")
  if (test.use %in% methods.noprefiliter) {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  # error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  # feature selection (based on percentages)
  data <- switch(
    EXPR = slot,
    'scale.data' = counts,
    object
  )
  if (is.null(x = reduction)) {
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(x = data.frame(data[features, cells.1, drop = FALSE]) > thresh.min) /
        length(x = cells.1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data.frame(data[features, cells.2, drop = FALSE]) > thresh.min) /
        length(x = cells.2),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > min.pct))
    if (length(x = features) == 0) {
      stop("No features pass min.pct threshold")
    }
    alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
    features <- names(
      x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
    )
    if (length(x = features) == 0) {
      stop("No features pass min.diff.pct threshold")
    }
  } else {
    data.alpha <- data.frame(
      pct.1 = rep(x = NA, times = length(x = features)),
      pct.2 = rep(x = NA, times = length(x = features))
    )
  }
  # feature selection (based on average difference)
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
      },
      function(x) {
        return(log(x = mean(x = x) + pseudocount.use))
      }
    )
  } else {
    mean
  }
  data.1 <- apply(
    X = data[features, cells.1, drop = FALSE],
    MARGIN = 1,
    FUN = mean.fxn
  )
  data.2 <- apply(
    X = data[features, cells.2, drop = FALSE],
    MARGIN = 1,
    FUN = mean.fxn
  )
  total.diff <- (data.1 - data.2)
  if (is.null(x = reduction) && slot != "scale.data") {
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff > logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) > logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      stop("No features pass logfc.threshold threshold")
    }
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    # Should be cells.1 and cells.2?
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }

  return(features)
}
