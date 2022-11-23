"
Create report of integration meticss

Usage:
    integration-metrics.R [options] <file>
    
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
  out_dir <- stringr::str_replace(
    stringr::str_replace(in_file, "data", "analysis"), ".h5ad", "/"
  )
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_file, layers = FALSE)
  
  # Plot label-free metrics ----------------------------------------------------
  data <- read_slot_h5ad(in_file, "labelfree")
  data$method <- row.names(data)
  data$method[data$method == "hlca"] <- c("counts_2000_scArches_HLCA")
  
  # Remove rows
  data <- data[which(!data$method == "X_emb"), ]
  
  df <- tidyr::gather(data, "metric", "score", -method)
  df$method <- factor(df$method, data$method)
  
  df$input <- stringr::str_c(
    stringr::str_split(df$method, "_", simplify = TRUE, n = 3)[, 1], "_",
    stringr::str_split(df$method, "_", simplify = TRUE, n = 3)[, 2]
  )
  df$algorithm <- stringr::str_split(df$method, "_", simplify = TRUE, n = 3)[, 3]
  
  ggplot2::ggplot(df, ggplot2::aes(algorithm, score)) +
    ggplot2::geom_point(ggplot2::aes(col = algorithm, shape = input), 
                        size = 5, position = "jitter") +
    ggplot2::facet_wrap(~metric, scales = "free") +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_classic(20) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_blank(), legend.margin = ggplot2::margin()
    )
  
  fn <- paste0(out_dir, "integration_metrics.png")
  ggplot2::ggsave(fn, height = 5.5, width = 12)
  
  # Label-metrics --------------------------------------------------------------
  cols <- names(ds@metadata$metrics)
  cols <- cols[cols != "labelfree"]
  data <- list()
  for (i in cols) {
    data[[i]] <- as.data.frame(ds@metadata$metrics[[i]])
    data[[i]]$type <- i
  }
  data <- dplyr::bind_rows(data)
  data$method <- data$X_index
  data$X_index <- NULL
  data$method[data$method == "hlca"] <- c("counts_2000_scArches_HLCA")
  
  # Remove rows
  data <- data[which(!data$method %in% c("X_emb", "counts_SCoV2")), ]
  
  df <- tidyr::gather(data, "metric", "score", -method, -type)
  df$method <- factor(df$method, data$method)
  
  df$input <- stringr::str_c(
    stringr::str_split(df$method, "_", simplify = TRUE, n = 3)[, 1], "_",
    stringr::str_split(df$method, "_", simplify = TRUE, n = 3)[, 2]
  )
  df$algorithm <- stringr::str_split(df$method, "_", simplify = TRUE, n = 3)[, 3]
  df$label <- c("TRUE" = "Control", "FALSE" = NA)[as.character(
    stringr::str_detect(df$algorithm, "PCA"))]
  
  ggplot2::ggplot(df, ggplot2::aes(algorithm, score)) +
    ggplot2::geom_point(ggplot2::aes(col = algorithm, shape = input), 
                        size = 5, position = "jitter") +
    ggplot2::facet_wrap(metric~type, scales = "free", ncol = 4) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_classic(20) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_blank(), legend.margin = ggplot2::margin()
    )
  
  fn <- paste0(out_dir, "integration_metrics-celltype.png")
  ggplot2::ggsave(fn, height = 6, width = 14)
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}