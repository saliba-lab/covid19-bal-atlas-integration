#' Summarize groups of columns in a matrix
#'
#' @param x A matrix-like object
#' @param groups A factor of group labels equal to ncol(x)
#' @param FUN A function to summarize columns
#' (default value requires sparseMatrixStats)
summarize_groups <- function (
    x = NULL, groups = NULL, FUN = sparseMatrixStats::rowMeans2
) {
  
  stopifnot(
    any(class(x) %in% c("matrix", "dgCMatrix", "dgTMatrix")),
    is.factor(groups) & length(groups) == ncol(x)
  )
  
  rn <- rownames(x)
  index <- split(1:length(groups), groups)
  
  fun <- function(j) {
    if (length(j) > 1) {
      FUN(x[, j])
    } else {
      x[, j]
    }
  }
  matrix <- sapply(index, fun)
  
  rownames(matrix) <- rn
  
  return(matrix)
}

#' Cluster cells
#' 
#' @param ds Data set of class 'SingleCellExperiment
#' @param alg Algorithm for clustering ("louvain", "leiden")
#' @param emb Embedding in reducedDims(ds)
#' @param dims Number of dimensions to choose from embedding
#' @param res_... Parameters setting the resolution parameter
#' 
#' @returns Data.frame containing factors of group membership
#' 
cluster_cells <- function(ds = NULL, emb = "X_mnn", alg = "louvain",
                          dims = 30,
                          res_from = 0.1, res_to = 3, res_by = 0.1) {
  
  mat <- SingleCellExperiment::reducedDim(ds, emb)
  g <- scran::buildSNNGraph(t(mat[, 1:dims]))
  res <- seq(res_from, res_to, by = res_by)
  data <- data.frame(row.names = colnames(ds))
  alg <- "louvain"
  for (i in res) {
    key <- paste0(i)
    if (alg == "louvain") {
      data[[key]] <- factor(
        igraph::cluster_louvain(g, resolution = i)$membership
      )
    }
    if (alg == "leiden") {
      data[[key]] <- factor(
        igraph::cluster_leiden(g, resolution_parameter = i)$membership
      )
    }
  }
  
  return(data)
}

#' Calculate AUCell scores for cell type markers
#' 
#' @param ds Dataset of type 'SingleCellExperiment'
#' @param markers List of marker genes
#' 
#' @returns Data.frame
#' 
run_AUCell <- function(lds = NULL, markers = NULL, batch = "patient", 
                       assay = "X") {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    class(markers) == "list",
    !is.null(batch),
    batch %in% names(ds@colData)
  )
  
  data <- list()
  for (i in unique(ds[[batch]])) {
    # Get data
    index <- which(ds[[batch]] == i)
    mat <- ds@assays@data[[assay]][, index]
    dimnames(mat) <- list(rownames(ds), colnames(ds)[index])
    
    # Run AUCell
    cells_rankings <- AUCell::AUCell_buildRankings(mat)
    cells_AUC <- AUCell::AUCell_calcAUC(markers, cells_rankings)
    data[[i]] <- as.data.frame(t(cells_AUC@assays@data$AUC))
  }
  data <- dplyr::bind_rows(data)
  data <- data[match(rownames(data), colnames(ds)), ]
  
  return(data)
}

#' Calculate mean expression of gene markers across clusters
#' 
#' @param ds Data set of type 'SingleCellExperiment'
#' @param markers List of marker genes by type
#' 
#' @returns Data.frame
#' 
calc_expression_scores <- function(ds=NULL, markers=NULL, scale=FALSE,
                                   assay = "X") {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    class(markers) == "list"
  )
  
  # Get expression data
  mat <- ds@assays@data[[assay]]
  dimnames(mat) <- dimnames(ds)
  
  data <- data.frame(
    row.names = colnames(mat)
  )
  for (i in names(markers)) {
    print(paste("Adding", i))
    index <- markers[[i]]
    if (any(!index %in% rownames(ds))) {
      miss <- index[!index %in% rownames(ds)]
      miss <- stringr::str_flatten(miss, collapse = " ")
      message(paste("Markers", miss, "are missing"))
    }
    index[index %in% rownames(ds)]
    index <- index[index %in% rownames(ds)]
    data[[i]] <- Matrix::colMeans(mat[index, ])
    
    if (scale) {
      data[[i]] <- scale(data[[i]])[, 1]
    }
    
  }
  
  return(data)
}

#' Annotate cells
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' @param groups Factor of cluster labels
#' @param AUC Data.frame of cell type scores
#' 
annotate_cells <- function(scores=NULL) {
  
  stopifnot(
    class(scores) == "data.frame"
  )
  
  cells <- row.names(scores)
  scores$cell <- cells
  scores <- tidyr::gather(scores, "type", "score", -cell)
  scores <- dplyr::mutate(
    dplyr::group_by(scores, cell), max = max(score)
  )
  scores <- dplyr::filter(
    dplyr::group_by(scores, cell), score == max(score)
  )
  dup_cells <- unique(scores$cell[which(duplicated(scores$cell))])
  scores <- scores[!duplicated(scores$cell), ]
  scores$type[scores$cell %in% dup_cells] <- "Unknown"
  
  scores <- scores$type[match(cells, scores$cell)]
  
  return(scores)
}


#' Annotate cells by cluster
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' @param groups Factor of cluster labels
#' @param AUC Data.frame of cell type scores
#' 
annotate_cells_by_cluster <- function(ds=NULL, groups=NULL, res=0.1, scores=NULL,
                           nrow=2, pt.size=2,
                           pl.write=FALSE, pl.height=6, pl.width=12) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    groups %in% names(ds@metadata),
    scores %in% names(ds@metadata)
  )
  
  # Summarize AUC scores across groups
  clust <- ds@metadata[[groups]]
  clust <- clust[, stringr::str_detect(names(clust), as.character(res))]
  data <- summarize_groups(x = t(ds@metadata[[scores]]), groups = clust)
  
  # Plot scores across groups
  df <- as.data.frame(t(data))
  df$cluster <- factor(rownames(df), rownames(df))
  df <- tidyr::gather(df, "type", "score", -cluster)
  df <- dplyr::summarise(
    dplyr::group_by(df, cluster), type = type, score = score, 
    mean = mean(score)
  )
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(type, score, color = type)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = mean)) +
    ggplot2::geom_point(position = "jitter", size = 3) +
    ggplot2::facet_wrap(~cluster, scales = "free", nrow = nrow) +
    ggplot2::theme_classic(15) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 5)
      )
    )
  print(p1)
  
  if (pl.write) {
    fn <- paste0(out_dir, scores, "_", groups, res, ".", "png")
    ggplot2::ggsave(fn, p1, height=pl.height, width=pl.width)
  }
  
  classify_max <- function(x) {which(x == max(x))}
  
  # Create annotations
  ann <- rownames(data)[apply(data, 2, classify_max)]
  names(ann) <- colnames(data)
  annotation <- factor(ann[clust], levels = rev(rownames(data)))
  
  return(annotation)
}


#' Plot cell-based scores on embedding
#' 
#' @param ds Data set of class 'SingleCellExperiment
#' @param scores Data.frame with score for each cell
#' @param embedding Matrix of coordinates for each cell
#' @param nrow Number of rows for facet_wrap layout
#' 
#' @returns plot
#' 
plot_cell_scores <- function(ds = NULL, scores = NULL, embedding = NULL,
                             nrow = 1, dims = c(1,2), 
                             pt.size = .1, pt.shape = 21, pt.stroke = .1,
                             pl.write=FALSE, pl.height=6, pl.width=12) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    scores %in% names(ds@metadata),
    embedding %in% SingleCellExperiment::reducedDimNames(ds)
  )
  
  emb <- SingleCellExperiment::reducedDim(ds, embedding)[, dims]
  colnames(emb) <- c("x", "y")
  sco <- ds@metadata[[scores]]
  data <- cbind(emb, sco)
  data <- tidyr::gather(data, "type", "score", -x, -y)
  data <- data[order(data$score), ]
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, col = score)) +
    ggplot2::geom_point(size = pt.size, shape = pt.shape, stroke = pt.stroke) +
    viridis::scale_color_viridis(option = "A", direction = -1) +
    ggplot2::facet_wrap(~type, nrow = nrow) +
    ggplot2::coord_fixed() +
    ggplot2::theme_light(base_size = 15) +
    ggplot2::labs(x = "UMAP-1", y = "UMAP-2") +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(barheight = 15, barwidth = 1, ticks = FALSE)
    )
  
  if (pl.write) {
    fn <- paste0(out_dir, embedding, "_", scores, ".", "png")
    ggplot2::ggsave(fn, plot, width=pl.width, height=pl.height)
  }
  
  return(plot)
}

#' Draw arrow for embedding
#' 
#' @param data Data.frame with x and y coordinates
#' @param direction Direction of the arrow (up, right)
#' 
#' @returns arrow
arrow <- function(data=NULL, direction="up") {
  
  stopifnot(
    !is.null(df)
  )
  
  right <-  ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = seq(min(data$x), max(data$x), length.out = 100)[25], 
    y     = min(data$y)*1.1, 
    yend  = min(data$y)*1.1, 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )
  
  up <- ggplot2::annotate(
    geom  = "segment", 
    x     = min(data$x)*1.1, 
    xend  = min(data$x)*1.1, 
    y     = min(data$y)*1.1, 
    yend  = seq(min(data$y), max(data$y), length.out = 100)[25], 
    arrow = grid::arrow(
      angle  = 30, length = ggplot2::unit(0.2, "inches"), type   = "closed"
    )
  )
  
  if (direction == "up") {return(up)}
  if (direction == "right") {return(right)}
}

#' Summarize rows that have overlapping coordinates
#' 
#' @param df Data.frame
#' @param x Column key
#' @param y Column key
#' @param FUN Function to summarize by
#' @param decimals Integer of decimals to round to
#' 
#' @returns Data.frame
#' 
summarize_overlapping_rows <- function(df=NULL, x="x", y="y", FUN=mean,
                                       breaks = 1000) {
  
  stopifnot(
    class(df) == "data.frame",
    class(df[[x]]) == "numeric",
    class(df[[y]]) == "numeric"
  )
  
  # Round coordinates
  for (i in c(x, y)) {
    df[[i]] <- as.numeric(cut(df[[i]], breaks))
  }
  df$match <- paste(df$x, df$y, sep = ":")
  
  # Summarize
  df <- dplyr::summarise(
    dplyr::group_by(df, match),
    x = mean(x), y = mean(y), col = FUN(col)
  )
  
  # Return as data.frame
  return(as.data.frame(df))
}


#' Plot embedding
#' 
#' @param ds Dataset of type 'SingleCellExperiment'
#' @param color Gene or coldata to color cells by
#' @param embedding Matrix with coordinates
#' @param dims Dimensions to plot from embedding
#' @param cells Barcodes of cells to plot
#' @param pt.size Point size
#' @param pt.aggr Boolean whether to aggregate overlapping points
#' @param pt.shape Point shape
#' @param pt.stroke Point stroke
#' @param pt.alpha Point alpha
#' 
#' @returns plot
#' 
plot_embedding <- function(ds=NULL, color=NA,
                           embedding=names(ds@int_colData$reducedDims)[1], 
                           assay = "X", cells = colnames(ds),
                           dims=1:2, pt.size=.1, pt.shape=20, 
                           pt.stroke=.1, pt.alpha=1, 
                           pl.write=FALSE, pl.height=6, pl.width=12
) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment"
  )
  
  df <- ds@int_colData$reducedDims[[embedding]][, dims]
  df <- as.data.frame(df)
  names(df) <- c("x", "y")
  
  if (color %in% rownames(ds)) {
    index <- match(color, rownames(ds))
    df$col <- ds@assays@data[[assay]][index, ]
  }
  if (color %in% names(ds@colData)) {
    df$col <- ds[[color]]
  }
  if (is.null(df$col)) {
    warning("Color key '", color, "' not found.")
    df$col <- NaN
  }
  
  if (class(df$col) %in% c("numeric", "integer", "array")) {
    colorscale <- viridis::scale_color_viridis(option = "A", direction = -1)
    guides <- ggplot2::guides(
      color = ggplot2::guide_colorbar(
        barheight = 20, barwidth = 1, ticks = FALSE
      )
    )
    pt.aggr <- TRUE
  }
  if (class(df$col) %in% c("factor", "character", "boolean", "logical")) {
    colorscale <- ggplot2::scale_color_discrete()
    guides <- ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 8, shape = 20)
      )
    )
    pt.aggr <- FALSE
  }
  
  # Subset
  index <- match(cells, colnames(ds))
  df <- df[index, ]
  
  if (pt.aggr) {
    # Set aggregation params
    if (all(is.nan(df$col))) {
      df$col <- 1
      FUN <- sum
      color <- "Density (cells/dot)"
    } else {
      FUN <- mean
    }
    # Summarize
    df <- summarize_overlapping_rows(df, breaks = 1000, FUN=FUN)
  }
    
  
  # Re-order
  df <- df[order(df$col), ]
  
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = col)) +
    ggplot2::geom_point(size=pt.size, shape=pt.shape, stroke=pt.stroke, 
                        alpha=pt.alpha) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = color, color = NULL) +
    ggplot2::theme_void(20) +
    colorscale + guides +
    arrow(df, "up") + arrow(df, "right")
  
  if(pl.write) {
    fn <- paste0(out_dir, embedding, "_", color, ".", "png")
    ggplot2::ggsave(fn, height=pl.height, width=pl.width, bg="white")
  }
  
  return(plot)
}

#' Plot cell-clusters on embedding
#' 
#' @param ds Data set of class 'SingleCellExperiment
#' @param clusters Key for cluster matrix in 'metadata'
#' @param embedding Key for reducedDims(ds)
#' @param nrow Number of rows for facet_wrap layout
#' 
#' @returns plot
#' 
plot_clusters_embedding <- function(ds=NULL, clusters=NULL, 
                                    embedding=NULL, dims=c(1,2),
                                    rt.plot=TRUE, write=FALSE, out_dir=NULL,
                                    pl.width = 20, pl.height = 9,
                                    nrow=3, ncol_legend=2, pt.size=.1, 
                                    pt.shape=21, pt.stroke=.1,
                                    lb.size=1.5, lb.pad=.1, lb.alpha=.4) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    clusters %in% names(ds@metadata),
    embedding %in% SingleCellExperiment::reducedDimNames(ds),
    !(write & is.null(out_dir))
  )
  
  # Fetch data
  emb <- SingleCellExperiment::reducedDims(ds)[[embedding]][, dims]
  colnames(emb) <- c("x", "y")
  clust <- ds@metadata[[clusters]]
  data <- cbind(emb, clust)
  data <- tidyr::gather(data, "type", "score", -x, -y)
  data$type <- factor(as.numeric(data$type))
  
  # Get label positions
  label_pos <- dplyr::group_by(data, type, score)
  label_pos <- dplyr::summarise(label_pos, x = mean(x), y = mean(y))
  
  # Plot
  plot <- ggplot2::ggplot(data, ggplot2::aes(x, y, col=score, label=score)) +
    ggplot2::geom_point(size=pt.size, shape=pt.shape, stroke=pt.stroke) +
    #ggrepel::geom_text_repel(data = label_pos, size = lb.size) +
    ggplot2::geom_label(data = label_pos, size = lb.size, 
                        label.padding = ggplot2::unit(lb.pad, "lines"),
                        alpha = lb.alpha) +
    ggplot2::facet_wrap(~type, nrow = nrow) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "grey90")
    ) +
    ggplot2::labs(x = "UMAP-1", y = "UMAP-2") +
    ggplot2::guides(
      color = ggplot2::guide_none()
    )
  
  if (write) {
    fn <- paste0(out_dir, embedding, "_", "cluster", "_", clusters, ".", "png")
    ggplot2::ggsave(fn, width = pl.width, height = pl.height)
  }
  if (rt.plot) {
    return(plot)
  }
  
}

#' Annotate clusters for each batch
#' 
#' @param sce SingleCellExperiment
#' @param batch Column in colData(sce)
#' @param assay Slot in assays(sce)
#' 
#' @returns Named vector with cell labels
#' 
per_sample_annotation <- function(ds=NULL, sample="patient", assay="cp10k",
                                  labels = NULL, resolution = 2,
                                  nv = 30, k = 10) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    sample %in% names(ds@colData),
    assay %in% names(ds@assays@data),
    class(labels) == "data.frame"
  )
  
  # Create object
  annotation <- list()
  
  # Add per-sample cell labels
  for (i in levels(ds[[sample]])) {
    print(paste("Running annotation for sample", i))
    set.seed(42)
    
    # Filter cells
    cells <- colnames(ds)[which(ds[[sample]] %in% i)]
    sce <- ds[, cells]
    sce@colData <- S4Vectors::DataFrame(labels[cells, ])
    
    # Filter HVGs
    hvg <- scran::modelGeneVar(sce, assay.type = assay)
    hvg <- scran::getTopHVGs(hvg, n = 2000)
    sce <- sce[hvg, ]
    
    # PCA
    sce <- scater::runPCA(sce, ncomponents = nv, exprs_values = "X")
    sce <- scater::runUMAP(sce, dimred="PCA")
    
    # KNN graph
    knn <- bluster::makeSNNGraph(sce@int_colData$reducedDims$PCA, k = k)
    
    # Cluster
    sce$cluster <- igraph::cluster_leiden(knn, resolution_parameter = resolution
                                          )$membership
    length(unique(sce$cluster))
    
    # Annotate
    ann <- as.data.frame(sce@colData)
    ann <- tidyr::gather(ann, "label", "value", -cluster)
    ann$N <- 1
    ann <- dplyr::summarize(
      dplyr::group_by(ann, cluster, label, value), N = sum(N)
    )
    ann <- dplyr::mutate(dplyr::group_by(ann, cluster, label), total = sum(N))
    ann$frac <- round(ann$N / ann$total, 2)
    
    ggplot2::ggplot(ann, ggplot2::aes(cluster, frac, fill=value)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::theme_classic(15) +
      ggplot2::facet_wrap(~label)
    
    ann <- dplyr::mutate(
      dplyr::group_by(ann, cluster, label), max = max(frac)
    )
    ann <- dplyr::filter(
      dplyr::group_by(ann, cluster, label), frac == max
    )
    
    # TODO: Break ties instead of discarding rows
    ann <- ann[!ann$cluster %in% which(
      matrixStats::rowMaxs(table(ann$cluster, ann$label)
                           ) > 1), ]
    
    # Assign labels
    ann <- ann[, c("cluster", "label", "value")]
    ann <- tidyr::spread(ann, label, value)
    ann <- ann[match(sce$cluster, ann$cluster), ]
    ann <- data.frame(ann)
    ann[is.na(ann)] <- "Unlabeled"
    rownames(ann) <- colnames(sce)
    
    # Add annotation
    ann$cluster <- NULL
    annotation[[i]] <- ann
  }
  
  # Combine annotations across sample
  annotation <- dplyr::bind_rows(annotation)[colnames(ds), ]
  
  return(annotation)
}