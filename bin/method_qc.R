"
Quality control

Usage:
    analysis_quality-control.R [options] <file>
    
Options:
    -h --help   Show this screen.
" -> doc

# Source functions
suppressMessages({
  for (i in list.files("bin/functions", full.names = TRUE)) {
    source(i)
  }
})

#' Add quality metrics
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' @param csv CSV file to add sample name based on barcode
#' @param xls Excel file with sample metadata
#' @param assay Name of main assay
#' 
#' @returns The function output
add_quality_metrics <- function(ds = NULL, assay = "X") {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    !is.null(libs),
    !is.null(samples),
    assay %in% names(ds@assays)
  )
  
  # Add QC metrics
  ds$libsize <- Matrix::colSums(ds@assays@data[[assay]])
  ds$nfeatures <- Matrix::colSums(ds@assays@data[[assay]] > 0)
  index <- which(stringr::str_detect(rownames(ds), "^MT-"))
  ds$percent.mt <- round(
    Matrix::colSums(ds@assays@data[[assay]][index, ]) / ds$libsize, 3
  ) * 100
  
  return(ds)
}

#' Assign quality loess
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' @param 
#' @param loess.t Numeric threshold for loess model
#' @param mad Numeric specifying the MAD used for adaptive thresholds
#' @param min_... Numeric lower bound for genes, counts, % mito. (mt)
#' @param max_... Numeric upper bound  for genes, counts, % mito. (mt)
#' 
#' @returns The function output
assign_quality_loess <- function(
    ds, loess.t = -0.15, span = 1
) {
  
  stopifnot(!is.null(ds))
  
  # Label conversion
  b2l <- c("TRUE" = "good", "FALSE" = "bad")
  
  # Loess model
  # Loess fit
  fit <- log10(as.data.frame(ds@colData[, c("libsize", "nfeatures")]))
  fit <- loess(nfeatures ~ libsize, fit, span = span)
  df <- as.data.frame(fit[c("x", "y", "fitted", "residuals")])
  
  # Model fit
  plot <- ggplot2::ggplot(df, ggplot2::aes(libsize, y)) +
    ggplot2::geom_point(shape = 1, alpha = .5, size = .5) +
    ggplot2::geom_line(ggplot2::aes(y = fitted), col = "blue")
  print(plot)
  
  # Loess threshold
  plot <- ggplot2::ggplot(df, ggplot2::aes(libsize, residuals)) +
    ggplot2::geom_point(shape = 1, alpha = .5, size = .5) +
    ggplot2::geom_line(ggplot2::aes(y = 0), col = "blue") +
    ggplot2::geom_line(ggplot2::aes(y = loess.t), col = "red")
  print(plot)
  
  # Store quality assignments
  label <- b2l[as.character(fit$residuals > loess.t)]
  
  return(label)
}

#' Assign quality by linear thresholds
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' @param min_... Numeric lower bound for genes, counts, % mito. (mt)
#' @param max_... Numeric upper bound  for genes, counts, % mito. (mt)
#' 
#' @returns Factor with quality assignment
#' 
assign_quality_linear <- function(
    ds,
    min_genes = 0, max_genes = Inf,
    min_counts = 0, max_counts = Inf, 
    min_mito = 0, max_mito = Inf
) {
  
  stopifnot(!is.null(ds))
  
  # Label conversion
  b2l <- c("TRUE" = "good", "FALSE" = "bad")
  
  # Linear thresholds
  label <- b2l[as.character(
    dplyr::between(ds$libsize, min_counts, max_counts) &
      dplyr::between(ds$nfeatures, min_genes, max_genes) & 
      dplyr::between(ds$percent.mt, min_mito, max_mito)
  )]
  
  return(label)
}

#' Assign quality
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' 
#' @returns Factor with quality labels
#' 
assign_quality_adaptive <- function(ds=NULL) {
  
  stopifnot(!is.null(ds))
  
  # Label conversion
  b2l <- c("TRUE" = "good", "FALSE" = "bad")
  
  # Adaptive thresholds
  df <- scuttle::perCellQCFilters(
    ds@colData, batch = ds$sample,
    sum.field = "libsize", detected.field = "nfeatures", 
    sub.fields = c("percent.mt")
  )
  label <- b2l[as.character(!df$discard)]
  
  return(label)
}


#' Violin plot of single QC metrics
#' 
#' @param ds Dataset of type 'SingleCellExperiment'
#' @param group Column key for grouping vector
#' @param metrics Column keys for QC metrics (numeric)
#' @param log_y Whether to log10-transform axes (boolean)
#' @param nrow Number of rows for faceting
#' 
#' @returns Plot
#' 
plot_qc_violin <- function(ds = NULL, color = NULL,
                           metrics = c("libsize", "nfeatures", "percent.mt"),
                           log_y = FALSE, nrow = 1, box.width = 0.1,
                           pl.write=FALSE, pl.width=8, pl.height=6, pl.dir=NULL
) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment"
  )
  
  # Log-transform axes
  log_y <- if(log_y) {ggplot2::scale_y_continuous(trans = "log10")} else {NULL}
  
  # Data Wrangling
  df <- as.data.frame(ds@colData[, c(metrics)])
  if (is.null(color)) {
    df$method <- ""
  } else {
    df$method <- factor(ds[[color]])
  }
  df <- tidyr::gather(df, "metric", "value",-method)
  
  # Plot
  plot <- ggplot2::ggplot(
    df, ggplot2::aes(x = method, y = value)
  ) + 
    ggplot2::geom_boxplot(width = box.width) +
    ggplot2::geom_violin(
      scale = "width", size = 1, fill = NA
    ) +
    ggplot2::facet_wrap(~metric, scales = "free", nrow = nrow) +
    log_y +
    ggplot2::theme_classic(20) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.y = ggplot2::element_line(size = 1),
      panel.grid.minor.y = ggplot2::element_line(size = .5)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(shape = 20, size = 5)
      )
    )
  
  if (pl.write) {
    fn <- paste0(pl.dir, "violin", "_", "metrics", "_", color, ".", "png")
    ggplot2::ggsave(fn, width = pl.width, height = pl.height)
  }
  
  return(plot)
}

#' Plot count-metric association by method
#'
#' @param ds Dataset of class 'SingleCellExperiment'
#' @param methods Column keys indicating methods to assign quality labels
#' @param sample_col Column key indicating samples
#' @param sample Sample(s) to plot (from sample_col)
#' @param x Numeric to compare metrics against
#' @param metrics Column keys for metrics to compare against x
#' @param cells Character of barcodes to plot
#' @param ncol Numeric number of columns for layout
#' @param log_... Boolean whether to log10-transform axes
#' 
#' @returns Plot
#' 
plot_count_association <- function(
    ds=NULL, method=NULL, sample_col="sample", sample=NULL, 
    x="libsize", metrics=c("nfeatures", "percent.mt"),
    log_x = TRUE, log_y = TRUE, 
    pl.dir=NULL, pl.write=FALSE, pl.height=6, pl.width=7
) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment",
    metrics %in% names(ds@colData)
  )
  
  # Log-transform axes
  log_x <- if (log_x) {ggplot2::scale_x_continuous(trans = "log10")} else {NULL}
  log_y <- if (log_y) {ggplot2::scale_y_continuous(trans = "log10")} else {NULL}
  
  # Fetch data
  df <- as.data.frame(
    ds@colData[, c(x, metrics)]
  )
  if (is.null(method)) {
    df$method <- "overall"
  } else {
    df$method <- ds[[method]]
  }
  
  # Subset by sample
  if (!is.null(sample) && sample %in% ds[[sample_col]]) {
    index <- ds[[sample_col]] == sample
    df <- df[index, ]
  }
  
  # Re-format
  df <- tidyr::gather(df, "metric", "value", -dplyr::all_of(x), -method)
  
  # Plot
  plot <- ggplot2::ggplot(df, ggplot2::aes(libsize, value)) +
    ggplot2::geom_bin_2d(bins = 100) +
    ggplot2::facet_grid(metric~method, scales = "free_y") +
    viridis::scale_fill_viridis(trans = "log10") +
    ggplot2::theme_classic(20) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.grid.major = ggplot2::element_line(size = 1),
      axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1)
    ) +
    ggplot2::labs(
      x = "Library size", y = NULL, fill = "cell/bin", title = sample
      ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        barheight = 15, barwidth = 1, ticks = FALSE
      )
    ) +
    log_x + log_y
  
  if (pl.write) {
    fn <- paste0(pl.dir, "scatter", "_", "metric", "_", method, ".", "png")
    ggplot2::ggsave(fn, width = pl.width, height = pl.height)
  }
  
  return(plot)
}


#' Plot association by sample
#' 
#' @param data Data.frame containing the cell data
#' @param color Column key to color cells
#' @param sample Column key with sample IDs
#' @param x Column key for x axis
#' @param y Column key for y axis
#' @param cells Character of barcodes to plot
#' @param log Boolean whether to log-transform axes
#' 
#' @returns Plot
#' 
plot_association_by_sample <- function(
    ds=NULL, quality=NULL, sample="sample", 
    x = "libsize", y = "nfeatures",
    log_x = TRUE, log_y = TRUE, f.row=4,
    pl.dir=NULL, pl.write=FALSE, pl.height=6, pl.width=12
) {
  
  stopifnot(
    class(ds) == "SingleCellExperiment"
  )
  
  # Log-transform axes
  log_x <- if (log_x) {ggplot2::scale_x_continuous(trans = "log10")} else {NULL}
  log_y <- if (log_y) {ggplot2::scale_y_continuous(trans = "log10")} else {NULL}
  
  # Subset data
  df <- as.data.frame(
    ds@colData[, c(x, y, sample)]
  )
  names(df) <- c("x", "y", "sample")
  if (is.null(quality)) {
    df$quality <- "overall"
  } else {
    df$quality <- ds[[quality]]
  }
  
  # Plot
  plot <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_bin2d(bins = 50) +
    viridis::scale_fill_viridis(trans = "log10") +
    ggplot2::facet_wrap(sample~quality, nrow = f.row) +
    ggplot2::labs(x = x, y = y) +
    ggplot2::theme_classic(15) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1)
    ) +
    log_x + log_y
  
  if (pl.write) {
    fn <- paste0(pl.dir, "scatter", "_", x, "_", y, "_", method, ".", "png")
    ggplot2::ggsave(fn, width = pl.width, height = pl.height)
  }
  
  return(plot)
}

#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  args <- docopt::docopt(doc)
  
  in_file <- args[["<file>"]]
  in_file_alt <- paste0(dirname(in_file), "/alt_", basename(in_file))
  in_file_xls <- "docs/overview.xlsx"
  h5file <- "data/BCB/filtered_feature_bc_matrix.h5"
  
  out_file <- stringr::str_replace(in_file, ".h5ad", "_qc_colData.csv")
  
  out_dir <- stringr::str_replace(in_file, "data", "analysis")
  out_dir <- stringr::str_replace(out_dir, "raw.h5ad", "qc/")
  
  if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}
  if (!endsWith(out_dir, "/")) {out_dir <- paste0(out_dir, "/")}
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_file)
  vs <- read_h5ad(in_file_alt)
  samples <- readxl::read_excel(in_file_xls, "samples")
  libs <- readxl::read_excel(in_file_xls, "libraries")
  
  # Add library info -----------------------------------------------------------
  
  ds@colData <- cbind(ds@colData, libs[match(ds$sample_id ,libs$libname), ])
  
  # Add virus counts
  index <- rownames(vs)
  index <- index[stringr::str_detect(index, "Hashtag", negate = TRUE)]
  for (i in index) {
    key <- paste0("virus_", i)
    ds[[key]] <- as.numeric(vs[i, ]@assays@data$X)
  }
  
  # Add hashtag counts
  index <- rownames(vs)
  index <- index[stringr::str_detect(index, "Hashtag", negate = FALSE)]
  ht <- t(as.matrix(vs[index, ]@assays@data$X))
  dimnames(ht) <- list(colnames(ds), index)
  
  
  # Demultiplex samples from libraries -----------------------------------------
  
  ind <- libs$libname[which(!is.na(libs$`TotalSeq-A`))]
  for (i in ind) {
    index <- libs$sample[libs$libname == i]
    index <- stringr::str_split(index, ",", simplify = TRUE)[1, ]
    
    index <- stringr::str_split(index, ":", simplify = TRUE)
    index <- as.data.frame(index)
    colnames(index) <- c("Hashtag", "Sample")
    index$Hashtag <- paste0("Hashtag", index$Hashtag, "_TotalA")
    
    cells <- colnames(ds)[ds$libname == i]
    df <- ht[cells, index$Hashtag]
    
    # TODO: Make function to assign hashtag labels
    
  }
  
  # Add sample info ------------------------------------------------------------
  
  # Add metrics ----------------------------------------------------------------
  message("Adding metadata & quality metrics")
  
  ds <- add_quality_metrics(ds, assay = "X")
  
  # Assign quality labels
  qc <- data.frame(
    row.names = rownames(ds@colData)
  )
  qc$linear <- assign_quality_linear(ds, min_genes = 300, max_mito = 15)
  qc$adaptive <- assign_quality_adaptive(ds)
  
  # Add columns to ds
  for (i in names(qc)) {
    key <- paste0("qc_", i)
    ds[[key]] <- factor(qc[[i]])
  }
  label_assign <- paste0("qc_", names(qc))
  
  
  # Plotting -------------------------------------------------------------------
  
  plot_qc_violin(ds, log_y = TRUE, pl.dir = out_dir, pl.write = TRUE)
  # TODO: Add filter option (e.g. qc_adaptive to look at plot after filtering)
  
  plot_qc_violin(ds, color = "sample_id", log_y = TRUE, nrow = 3,
                 pl.dir = out_dir, pl.write = TRUE, pl.height = 12, pl.width=12)
  plot_qc_violin(ds, color = "sample_id", log_y = TRUE, nrow = 3,
                 pl.dir = out_dir, pl.write = TRUE, pl.height = 8)
  
  plot_count_association(ds, pl.dir=out_dir, pl.write=T, pl.width=6)
  for (i in label_assign) {
    plot_count_association(ds, method = i,
                                     pl.dir=out_dir, pl.write=T, pl.width=7)
  }
  
  key <- "sample_id"
  ds[[key]] <- factor(ds[[key]])
  pdf(file = paste0(out_dir, "scatter_metric_sample.pdf"))
  for (i in levels(ds[[key]])) {
    print(plot_count_association(ds, sample = i, sample_col = key))
  }
  dev.off()
  
  for (i in label_assign) {
    pdf(file = paste0(out_dir, "scatter_metric_", i, "_", key,".pdf"))
    for (j in levels(ds[[key]])) {
      print(plot_count_association(ds, method=i, sample=j, sample_col=key))
    }
    dev.off()
  }
  
  # Write output ---------------------------------------------------------------
  message("Writing output to '", out_file, "'...")
  write.csv(as.data.frame(ds@colData), out_file)
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}