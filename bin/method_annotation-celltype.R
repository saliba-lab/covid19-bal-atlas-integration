"
Cell type annotation

Usage:
    annotation-markers.R --plot-dir=<path> [options] <file>
    
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
source("reports/sample-overview.R")

#' The main script function
main <- function() {
  
  # Libraries ------------------------------------------------------------------
  
  library('clustree')
  
  # Variables ------------------------------------------------------------------
  args <- docopt::docopt(doc)
  
  in_file <- args[["<file>"]]
  in_overview <- "docs/overview.xlsx"
  in_doublet <- stringr::str_replace(in_file, ".h5ad", "_doublets.tsv")
  in_hlca_labels <- stringr::str_replace(in_file, ".h5ad", "_hlca-labels.csv")
  in_hlca_latent <- stringr::str_replace(in_file, ".h5ad", "_hlca-latent.csv")
  
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
  ds <- read_h5ad(in_file, layers = "cp10k")
  
  # HLCA transfer
  hlca <- list()
  hlca$labels <- read.csv(in_hlca_labels, row.names = "X")
  
  # Overview
  ov <- list()
  for (i in readxl::excel_sheets(in_overview)) {
    ov[[i]] <- readxl::read_excel(in_overview, i)
  }
  ov <- clean_dataset_overview(ov)
  
  # Marker genes
  xls_url <- "https://nubes.helmholtz-berlin.de/s/id9rrKNeditKMmF"
  xls_url <- paste0(xls_url, "/", "download")
  xls <- "docs/celltype/markers.xlsx"
  download.file(xls_url, xls)
  markers <- list()
  for (i in readxl::excel_sheets(xls)) {
    markers[[i]] <- readxl::read_excel(xls, sheet = i)
  }
  
  # Public markers
  url <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
  panglao <- "docs/celltype/panglaoDB.csv"
  download.file(url, panglao)
  pan <- read.table(panglao, sep = "\t", header = TRUE)
  pan <- pan[stringr::str_detect(pan$species, "Hs"), ]
  
  # Published cohort 1
  url <- "https://ars.els-cdn.com/content/image/1-s2.0-S0092867421013830-mmc2.xlsx"
  cohort1 <- "docs/celltype/cohort1_macro.xlsx"
  download.file(url, cohort1)
  readxl::excel_sheets(cohort1)
  ch1 <- readxl::read_excel(cohort1, "2", skip = 2)
  ch1 <- ch1[ch1$FDR < 1e-15, ]
  
  # Create random vector to split dataset
  ds$random <- sample(1:10, ncol(ds), replace = TRUE)
  
  # Exploratory analysis -------------------------------------------------------
  
  # Select embedding
  names(ds@int_colData$reducedDims)
  emb <- "counts_seurat2000_scVI_n30l1h128"
  
  # Calculate UMAP embedding
  set.seed(42)
  key <- paste0(emb, "_", "umap")
  ds <- scater::runUMAP(ds, name=key, dimred=emb, metric = "cosine",
                        min_dist = 0.3)
  
  # Plot
  plot_embedding(ds, embedding = key, color = "")
  fn <- paste0(out_dir, "umap", "_", "density", ".", "png")
  ggplot2::ggsave(fn, width = 7, height = 6, bg = "white")
  
  # Plot
  plot_embedding(ds, embedding = key, color = "MKI67")
  fn <- paste0(out_dir, "umap", "_", "MKI67", ".", "png")
  ggplot2::ggsave(fn, width = 7, height = 6, bg = "white")
  
  # Plot
  plot_embedding(ds, embedding = key, color = "percent.mt")
  fn <- paste0(out_dir, "umap", "_", "percentMT", ".", "png")
  ggplot2::ggsave(fn, width = 7, height = 6, bg = "white")
  
  # Plot
  plot_embedding(ds, embedding = key, color = "libsize", pt.stroke = .5)
  fn <- paste0(out_dir, "umap", "_", "libsize", ".", "png")
  ggplot2::ggsave(fn, width = 7, height = 6, bg = "white")
  
  # Plot
  genes <- rownames(ds)[grep("IGH", rownames(ds))]
  fn <- paste0(out_dir, "IGH-genes.pdf")
  pdf(fn)
  for (i in genes) {
    p <- plot_embedding(ds, embedding = key, color = i)
    print(p)
  }
  dev.off()
  
  # Assign sample names to multiplexed libraries -------------------------------
  
  ht2lab <- c(
    "L36:1" = "doublet",
    "L36:2" = "BAL_44",
    "L36:3" = "BAL_45",
    "L36:4" = "BAL_45",
    "L37:1" = "BAL_47",
    "L37:2" = "BAL_47",
    "L37:3" = "BAL_46",
    "L37:4" = "BAL_46",
    "L39:1" = "BAL_32",
    "L39:2" = "negative",
    "L39:3" = "doublet",
    "L39:4" = "BAL_31",
    "L39:5" = "BAL_33",
    "L40:1" = "BAL_50",
    "L40:2" = "doublet",
    "L40:3" = "BAL_49",
    "L40:4" = "BAL_49",
    "L42:1" = "BAL_41",
    "L42:2" = "BAL_40",
    "L42:3" = "BAL_37",
    "L42:4" = "doublet",
    "L42:5" = "BAL_39",
    "L42:6" = "doublet"
  )
  
  index <- which(ds$hto_clust != "None")
  ds$sample <- as.character(ds$sample)
  ds$sample[index] <- ht2lab[ds$hto_clust[index]]
  lvl <- unique(ds$sample)[order(as.numeric(
    stringr::str_split(unique(ds$sample), "_", simplify = TRUE)[, 2]
    ))
  ]
  ds$sample <- factor(ds$sample, lvl)
  
  # Plot
  plot_embedding(ds, embedding = key, color = "sample", pt.stroke = .2)
  fn <- paste0(out_dir, "umap", "_", "sample", ".", "png")
  ggplot2::ggsave(fn, width = 10, height = 6, bg="white")
  
  # Assign patient data to samples ---------------------------------------------
  
  # Days post symptom onset
  ds$dpso <- ov$samples$dpso[match(ds$sample, ov$samples$sample)]
  
  # Plot
  plot_embedding(ds, embedding = key, color = "dpso")
  fn <- paste0(out_dir, "umap", "_", "dpso", ".", "png")
  ggplot2::ggsave(fn, width = 10, height = 6, bg="white")
  
  # Assign patient id
  s2p <- ov$samples$patient
  names(s2p) <- ov$samples$sample
  ds$patient <- factor(s2p[ds$sample])
  
  # Plot
  plot_embedding(ds, embedding = key, color = "patient", pt.stroke = .1)
  fn <- paste0(out_dir, "umap", "_", "patient", ".", "png")
  ggplot2::ggsave(fn, width = 10, height = 6, bg="white")
  
  # Patient outcome ------------------------------------------------------------
  
  ds$outcome <- "absent"
  ds$outcome[which(
    ds$patient %in% c("C19-CB-0916", "C19-CB-40001", "C19-CB-40003", "C19-CB-0915")
    )] <- "deceased"
  ds$outcome[which(
    ds$patient %in% c("C19-CB-0914", "C19-CB-40002", "C19-CB-0917", "C19-CB-0913")
  )] <- "survived"
  
  plot_embedding(ds, embedding = key, color = "outcome", pt.stroke = .1, 
                 pl.title = "Disease outcome", color.pal = "manual",
                 color.colors = c("grey", "indianred", "skyblue"))
  fn <- paste0(out_dir, "umap", "_", "outcome", ".", "png")
  ggplot2::ggsave(fn, width = 7.5, height = 6, bg="white")
  
  # Cell cycle annotation ------------------------------------------------------
  message("Annotating cell cycle...")
  
  # Compute AUC
  fn <- paste0(out_dir, "auc_cellcycle.csv")
  if (file.exists(fn)) {
    ds@metadata$cellcylce <- read.csv(fn, row.names = "X")
  } else {
    mrk <- markers$cellcycle
    mrk <- split(mrk$gene, mrk$type)
    ds@metadata$cellcycle <- run_AUCell(ds, mrk, batch = "random", assay = "X")
    names(ds@metadata$cellcycle) <- names(mrk)
    write.csv(ds@metadata$cellcycle, fn) 
  }
  
  # Plot
  plot_cell_scores(ds, scores = "cellcycle", embedding = key, nrow = 2)
  fn <- paste0(out_dir, "umap", "_", "cellcycle-scores", ".", "png")
  ggplot2::ggsave(fn, width = 4.2, height = 6)
  
  # Label
  data <- ds@metadata$cellcycle
  data$clust <- "G1"
  data$clust[data$G2M > 0.1 & data$S < 0.1] <- "G2M"
  data$clust[data$S > 0.1 & data$G2M < 0.1] <- "S"
  data$clust[data$S > 0.1 & data$G2M > 0.1] <- "S+G2M"
  ggplot2::ggplot(data, ggplot2::aes(G2M, S, col = clust)) +
    ggplot2::geom_point(size = .1, stroke = 1, shape = 16) +
    ggplot2::theme_classic(15) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5, shape = 20))
    )
  ds$cc1 <- factor(data$clust, c("G1", "S", "G2M", "S+G2M"))
  ds$cc2 <- "G1"
  ds$cc2[data$clust %in% c("S", "G2M", "S+G2M")] <- "Prolif."
  
  # Plot
  plot_embedding(ds, embedding = key, color = "cc1", pt.stroke = .1,
                 pl.title = "Cell cycle phase", color.pal = "manual", 
                 color.colors = c("grey", "navy", "indianred", "orange"))
  fn <- paste0(out_dir, "umap", "_", "cellcycle", ".", "png")
  ggplot2::ggsave(fn, width = 7.3, height = 6, bg = "white")
  
  # Clustering -----------------------------------------------------------------
  
  # Compute clusters
  fn <- paste0(
    out_dir, "cluster", "_", stringr::str_replace_all(emb,"_","-"), ".", "csv"
    )
  if (file.exists(fn)) {
    ds@metadata$cluster <- read.csv(fn, row.names = "X")
    names(ds@metadata$cluster) <- stringr::str_remove(
      names(ds@metadata$cluster), "X"
      )
  } else {
    ds@metadata$cluster <- cluster_cells(ds, emb = emb, alg = "leiden")
    write.csv(ds@metadata$cluster, fn)
  }
  
  # Create cluster tree
  data <- ds@metadata$cluster
  names(data) <- paste0("res_", names(data))
  
  tree <- clustree::clustree(data, prefix = "res_")
  tree
  fn <- paste0(out_dir, "clustree", ".", "png")
  ggplot2::ggsave(fn, width = 12, height = 12)
  
  # Plot
  ds$cluster <- factor(as.numeric(ds@metadata$cluster$`2.4`))
  plot_embedding(ds, embedding = key, color = "cluster")
  
  # Plot
  plot_clusters_embedding(ds, clusters="cluster", embedding=key)
  fn <- paste0(out_dir, "cluster", "-", "overview", ".", "png")
  ggplot2::ggsave(fn, width = 27, height = 9)
  
  # Cell type markers - level 1 ------------------------------------------------
  
  # Compute AUC
  fn <- paste0(out_dir, "auc_level1.csv")
  if (file.exists(fn)) {
    ds@metadata$level_1 <- read.csv(fn, row.names = "X")
  } else {
    mrk <- markers$level_1
    mrk <- split(mrk$gene, mrk$type)
    ds@metadata$level_1 <- run_AUCell(ds, mrk, batch = "random", assay = "X")
    names(ds@metadata$level_1) <- names(mrk)
    write.csv(ds@metadata$level_1, fn)
  }
  
  # Plot
  plot_cell_scores(ds, scores = "level_1", embedding = key, nrow = 2)
  fn <- paste0(out_dir, "umap", "_", "level-1-scores", ".", "png")
  ggplot2::ggsave(fn, width = 12, height = 6)
  
  # Cell type markers - level 2 ------------------------------------------------
  
  # Compute AUC
  fn <- "data/BCB/full_auc_level2.csv"
  if (file.exists(fn)) {
    ds@metadata$level_2 <- read.csv(fn, row.names = "X")
  } else {
    mrk <- markers$level_2
    mrk <- split(mrk$gene, mrk$type)
    ds@metadata$level_2 <- run_AUCell(ds, mrk, batch = "random", assay = "X")
    names(ds@metadata$level_2) <- names(mrk)
    write.csv(ds@metadata$level_2, fn) 
  }
  
  # Plot
  plot_cell_scores(ds, scores = "level_2", embedding = key, nrow = 2)
  fn <- paste0(out_dir, "umap", "_", "level-2-scores", ".", "png")
  ggplot2::ggsave(fn, width = 22, height = 6)
  
  # Cell type markers - level 3 ------------------------------------------------
  
  # Compute AUC
  fn <- "data/BCB/full_auc_level3.csv"
  if (file.exists(fn)) {
    ds@metadata$level_3 <- read.csv(fn, row.names = "X")
  } else {
    mrk <- markers$level_3
    mrk <- split(mrk$gene, mrk$type)
    ds@metadata$level_3 <- run_AUCell(ds, mrk, batch = "random", assay = "X")
    names(ds@metadata$level_3) <- names(mrk)
    write.csv(ds@metadata$level_3, fn) 
  }
  
  # Plot AUC
  plot_cell_scores(ds, scores = "level_3", embedding = key, nrow = 4)
  fn <- paste0(out_dir, "umap", "_", "level-3-scores", ".", "png")
  ggplot2::ggsave(fn, width = 16, height = 12, dpi = 100)
  
  # Cell type markers - PanglaoDB ----------------------------------------------
  
  # Compute AUC
  fn <- paste0(out_dir, "auc_panglao.csv")
  if (file.exists(fn)) {
    ds@metadata$panglao <- read.csv(fn, row.names = "X")
  } else {
    mrk <- split(pan$official.gene.symbol, pan$cell.type)
    ds@metadata$panglao <- run_AUCell(ds, mrk, batch = "random", assay = "X")
    names(ds@metadata$panglao) <- names(mrk)
    write.csv(ds@metadata$panglao, fn)
  }
  
  ds$score <- ds@metadata$panglao$Pulmonary.alveolar.type.I.cells
  plot_embedding(ds, embedding = key, color = "score")
  
  # Macrophage markers - cohort 1 ----------------------------------------------
  
  # Compute AUC
  fn <- "data/BCB/full_auc_cohort1-mp.csv"
  if (file.exists(fn)) {
    ds@metadata$cohort1_mp <- read.csv(fn, row.names = "X")
  } else {
    mrk <- split(ch1$gene, ch1$cluster)
    ds@metadata$cohort1_mp <- run_AUCell(ds, mrk, batch = "random", assay = "X")
    names(ds@metadata$cohort1_mp) <- names(mrk)
    write.csv(ds@metadata$cohort1_mp, fn)
  }
  
  # Plot
  plot_cell_scores(ds, scores = "cohort1_mp", embedding = key, nrow = 2)
  fn <- paste0(out_dir, "umap", "_", "cohort1-mp-scores", ".", "png")
  ggplot2::ggsave(fn, width = 10, height = 6)
  
  # Annotation - final ---------------------------------------------------------
  
  # Label cells
  data <- ds@metadata$level_2
  df <- data.frame(row.names = rownames(data))
  df$total <- rowSums(data)
  df$pos <- rowSums(data > attr(data, "thresh"))
  df$max <- sparseMatrixStats::rowMaxs(as.matrix(data))
  df$sd <- sparseMatrixStats::rowSds(as.matrix(data))
  
  index <- df$pos > 0
  data[index, ] == df$max[index]
  
  # Plot
  ds$score <- df$sd
  plot_embedding(ds, embedding = key, color = "score")
  
  
  df <- as.data.frame(auc)
  df$bc <- rownames(df)
  df <- tidyr::gather(df, "type", "score", -bc)
  df <- df[df$score > 0, ]
  ggplot2::ggplot(df, ggplot2::aes(type, score, col = type)) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::geom_point(
      data = data.frame(name = names(auc), thresh = attr(auc, "thresh")), 
      ggplot2::aes(name, thresh), color = "black"
      ) +
    ggplot2::theme_classic(15) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey70"),
      panel.grid.minor.y = ggplot2::element_line(color = "grey80"),
      axis.title.x = ggplot2::element_blank()
    )
  
  df$type <- stringr::str_replace(df$type, "\\.", "\\/")
  df$th <- attr(auc, "thresh")[df$type]
  df$keep <- df$score > df$th
  
  df <- df[df$keep, ]
  df <- dplyr::group_by(df, bc, type)
  df <- dplyr::summarise(df, score = score)
  
  # Explore
  names(auc)
  type <- names(auc)[3]
  hist(auc[[type]][auc[[type]] > 0], breaks = 100, main = type)
  ds$score <- auc[[type]] > 0.4
  plot_embedding(ds, embedding = key, cells = 30000, color = "score")
  
  plot_embedding(ds, embedding = key, cells = 30000, color = "JCHAIN")
  
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
