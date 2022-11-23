"
Cell type annotation

Usage:
    analysis_annotation-celltype.R --plot-dir=<path> [options] <file>
    
Options:
    -h --help             Show this screen.
    -o --out-dir=<path>   Output directory
" -> doc

# Source functions
suppressMessages({
  for (i in list.files("bin/functions", full.names = TRUE)) {
    source(i)
  }
})

#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  args <- docopt::docopt(doc)
  
  in_file <- args[["<file>"]]
  out_file <- in_file
  out_file_ann <- stringr::str_replace(out_file, ".h5ad", "_annotation.tsv")
  out_dir <- stringr::str_replace(
    stringr::str_replace(in_file, "data", "analysis"), ".h5ad", "/"
    )
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_file)
  
  # Select embedding
  names(ds@int_colData$reducedDims)
  emb <- "counts_seurat5000_scVI_n30l1h128"
  
  # Calculate UMAP embedding ---------------------------------------------------
  set.seed(42)
  key <- paste0(emb, "_", "umap")
  ds <- scater::runUMAP(ds, name=key, dimred=emb, metric = "cosine",
                        min_dist = 0.3)
  
  plot_embedding(ds, embedding = key, assay = "cp10k", color = "MKI67")
  
  # Explore HLCA predictions ---------------------------------------------------
  
  # Determine uncertainty threshold
  thresh <- 0.2
  cols  <- names(ds@colData)
  anns <- stringr::str_remove(cols[stringr::str_detect(cols, "pred")], "_pred")
  preds <- as.data.frame(ds@colData[, paste0(anns, "_pred")])
  uncs <- as.data.frame(ds@colData[, paste0(anns, "_uncertainty")])
  
  df <- tidyr::gather(uncs, "x", "y")
  df$x <- stringr::str_remove(df$x, "_uncertainty")
  ggplot2::ggplot(df, ggplot2::aes(x,y, col = x)) +
    ggplot2::geom_violin(scale = "width", size = 1) +
    ggplot2::geom_hline(yintercept = thresh) +
    ggplot2::theme_classic(20) +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "HLCA prediction", y = "Uncertainty")
  
  # Replace labels with exceeding uncertainty
  for (i in names(preds)) {preds[[i]] <- as.character(preds[[i]])}
  preds[uncs > thresh] <- "Unknown"
  
  # Per-sample annotation ------------------------------------------------------
  
  # HLCA predictions
  result <- per_sample_annotation(ds, labels = preds, resolution = 2)
  for (i in names(result)) {
    j <- stringr::str_remove(i, "_pred")
    ds[[i]] <- preds[[i]]
    ds[[j]] <- result[[i]]
  }
  
  # Marker gene annotation
  xls <- "docs/cell-type-markers.xlsx"
  preds <- data.frame(row.names = colnames(ds))
  for (i in c("level_1", "level_2")) {
    j <- paste0("marker_", i, "_pred")
    markers <- readxl::read_excel(xls, i)
    markers <- split(markers$gene, markers$type)
    scores <- run_AUCell(ds, markers = markers)
    preds[[j]] <- annotate_cells(scores)
  }
  
  result <- per_sample_annotation(ds, labels = preds, resolution = 2)
  for (i in names(result)) {
    j <- stringr::str_remove(i, "_pred")
    ds[[i]] <- preds[[i]]
    ds[[j]] <- result[[i]]
  }
  
  # Write results
  lbs <- as.data.frame(ds@colData[, c("ann_level_3", "ann_level_4")])
  fn <- paste0(out_dir, "labels.csv")
  write.csv(lbs, fn)
  
  # Annotate final-embedding ---------------------------------------------------
  
  knn <- bluster::makeSNNGraph(ds@int_colData$reducedDims[[emb]])
  ds$cluster <- igraph::cluster_leiden(knn, resolution_parameter=1)$membership
  
  ann <- data.frame(
    row.names = colnames(ds),
    pred = as.character(ds$ann_finest_level_pred),
    clust = ds$cluster
  )
  ann$N <- 1
  ann <- dplyr::summarize(
    dplyr::group_by(ann, clust, pred), N = sum(N)
  )
  ann <- dplyr::mutate(dplyr::group_by(ann, clust), total = sum(N))
  ann$frac <- round(ann$N / ann$total, 2)
  
  ggplot2::ggplot(ann, ggplot2::aes(clust, frac, fill=pred)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::theme_classic(15)
  
  ann <- dplyr::mutate(
    dplyr::group_by(ann, clust), max = max(frac)
  )
  ann <- dplyr::filter(
    dplyr::group_by(ann, clust), frac == max
  )
  ann <- ann[which(table(ann$clust) <= 1), ]
  ds$celltype <- ann$pred[match(ds$cluster, ann$clust)]
  ds$celltype[is.na(ds$celltype)] <- "Unknown"
  
  plot_embedding(ds, embedding = key, color = "celltype")
  
  # Clustering -----------------------------------------------------------------
  ds@metadata[["mnn_louvain"]] <- cluster_cells(ds, "X_mnn", "louvain", 
                                                res_to = 4)
  ds@metadata[["scVI_louvain"]] <- cluster_cells(ds, "X_scVI", "louvain", 
                                                 res_to = 4)
  
  plot_clusters_embedding(ds, "mnn_louvain", "X_mnn_umap", rt.plot = FALSE,
                          write = TRUE, out_dir = out_dir, nrow = 4,
                          pl.width = 22, pl.height = 9)
  
  plot_clusters_embedding(ds, "scVI_louvain", "X_mnn_umap", rt.plot = FALSE,
                          write = TRUE, out_dir = out_dir, nrow = 4,
                          pl.width = 22, pl.height = 9)
  
  plot_clusters_embedding(ds, "mnn_louvain", "X_scVI_umap", rt.plot = FALSE,
                          write = TRUE, out_dir = out_dir, nrow = 4,
                          pl.width = 16, pl.height = 9)
  
  plot_clusters_embedding(ds, "scVI_louvain", "X_scVI_umap", rt.plot = FALSE,
                          write = TRUE, out_dir = out_dir, nrow = 4,
                          pl.width = 16, pl.height = 9)
  
  # Cell cycle annotation ------------------------------------------------------
  markers <- Seurat::cc.genes.updated.2019
  names(markers) <- c("S-Phase", "G2M-Phase")
  ds@metadata[["CellCycle"]] <- calc_expression_scores(ds, markers)
  
  plot_cell_scores(ds, "CellCycle", embedding = "X_mnn_umap", nrow = 2,
                   pt.size = .1, pt.shape = 20, pt.stroke = .1,
                   pl.write = TRUE, pl.width = 4.5)
  
  plot_cell_scores(ds, "CellCycle", embedding = "X_scVI_umap", nrow = 2,
                   pt.size = .1, pt.shape = 20, pt.stroke = .1,
                   pl.write = TRUE, pl.width = 4.5)
  
  data <- ds@metadata$CellCycle
  names(data) <- c("x", "y")
  data$prolif <- data$x >= 0.5 | data$y >= 0.75
  
  ggplot2::ggplot(data, ggplot2::aes(x, y, col = prolif)) +
    ggplot2::geom_point(size = 1, shape = 21, strok = 1) +
    ggplot2::theme_classic(15)
  
  ds$Proliferating <- data$prolif
  plot_embedding(ds, "Proliferating", "X_mnn_umap", pt.size=.1, pt.stroke=.1)
  
  # Cell type annotation - level 1 ---------------------------------------------
  markers <- readxl::read_excel("docs/cell-type-markers.xlsx", "level_1")
  markers <- split(markers$gene, markers$type)
  
  # ds@metadata[["level_1_type"]] <- run_AUCell(ds, markers)
  ds@metadata$level_1_type <- calc_expression_scores(ds, markers, scale=FALSE)
  
  plot_cell_scores(ds, "level_1_type", nrow = 2,
                   embedding = "X_mnn_umap", 
                   pt.size = .1, pt.shape = 20, pt.stroke = .1,
                   pl.write = TRUE, pl.width = 12)
  
  plot_cell_scores(ds, "level_1_type", nrow = 2,
                   embedding = "X_scVI_umap", 
                   pt.size = .1, pt.shape = 20, pt.stroke = .1,
                   pl.write = TRUE, pl.width = 12)
  
  # Annotate based on MNN-integration
  ds$celltype_mnn <- annotate_cells(ds, groups = "mnn_louvain", res=3.5, 
                                    scores="level_1_type", 
                                    nrow = 2, pl.write=TRUE, pl.width=20)
  
  plot_embedding(ds, "celltype_mnn", "X_mnn_umap", pt.size=.1, pt.stroke=.2,
                 pl.write=TRUE, pl.width=8)
  plot_embedding(ds, "celltype_mnn", "X_scVI_umap", pt.size=.1, pt.stroke=.2,
                 pl.write=TRUE, pl.width=8)
  
  # Annotate based on scVI-integration
  ds$celltype_scvi <- annotate_cells(ds, groups = "scVI_louvain", res=3.5, 
                                    scores="level_1_type", 
                                    nrow = 2, pl.write=TRUE, pl.width=20)
  
  plot_embedding(ds, "celltype_scvi", "X_mnn_umap", pt.size=.1, pt.stroke=.2,
                 pl.write=TRUE, pl.width=8)
  plot_embedding(ds, "celltype_scvi", "X_scVI_umap", pt.size=.1, pt.stroke=.1,
                 pl.write=TRUE, pl.width=8)
  
  # Cell type annotation - ref_A -----------------------------------------------
  markers <- readxl::read_excel("docs/cell-type-markers.xlsx", "ref_A")
  m <- list()
  for (i in markers$type) {
    m[[i]] <- stringr::str_split(
      markers$genes[markers$type == i], ", "
      )[[1]]
  }
  
  # TODO: Find way to replace AUCell for large datasets
  # ds@metadata[["ref_A"]] <- run_AUCell(ds, m)
  ds@metadata[["ref_A"]] <- calc_expression_scores(ds, m)
  
  plot_cell_scores(ds, scores = "ref_A", embedding = "X_mnn_umap", nrow = 4,
                   pt.size = 1, pl.write = TRUE, pl.height = 9, pl.width = 18)
  
  plot_cell_scores(ds, scores = "ref_A", embedding = "X_scVI_umap", nrow = 4,
                   pt.size = 1, pl.write = TRUE, pl.height = 9, pl.width = 14)
  
  # Write output ---------------------------------------------------------------
  
  # View file
  file <- rhdf5::h5ls(out_file)
  file[file$group == "/obsm", ]
  
  # Add selected embeddings
  for (i in c("X_mnn_umap", "X_scVI_umap")) {
    print(paste("Adding", i))
    key <- paste0("obsm/", i)
    
    if (i %in% file$name) {
      rhdf5::h5delete(out_file, key)
    }
    
    obj <- as.matrix(t(ds@int_colData$reducedDims[[i]]))
    rhdf5::h5write(obj = obj, file = out_file, name = key)
  }
  
  # Add annotations
  obj <- data.frame(
    row.names = rownames(ds@colData)
  )
  for (i in c("celltype_mnn", "celltype_scvi", "Proliferating")) {
    obj[[i]] <- ds[[i]]
  }
  write.table(obj, file = out_file_ann, sep = "\t")
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}