"
Herpesvirus in COVID-19 BALs

Usage:
    analysis_annotation-celltype.R [options]
    
Options:
    -h --help             Show this screen.
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
  
  in_file <- "data/BCB/alt_raw.h5ad"
  in_full <- "data/BCB/full.h5ad"
  in_qc <- "data/BCB/raw_qc_colData.csv"
  
  out_dir <- "analysis/BCB/overview/"
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  
  ds <- read_h5ad(in_file)
  qc <- read.csv(in_qc, row.names = "X")
  
  # Add QC to object
  ds@colData <- cbind(ds@colData, qc)
  
  # Viral counts per sample ----------------------------------------------------
  message("Producing viral count summary...")
  
  # Select data
  cols <- names(ds@colData)
  cols <- cols[stringr::str_detect(cols, "virus")]
  data <- as.data.frame(ds@colData[, cols])
  names(data) <- stringr::str_remove(names(data), "virus_")
  data$lib <- ds@colData$sample
  data <- tidyr::gather(data, "virus", "count", -lib)
  
  # Summarize
  df <- dplyr::summarise(
    dplyr::group_by(data, lib, virus), count = sum(count)
    )
  
  # Plot
  ggplot2::ggplot(data, ggplot2::aes(virus, count, col = lib)) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::theme_light(15)
  fn <- paste0(out_dir, "viral-counts-violin", ".", "png")
  ggplot2::ggsave(fn, width = 12, height = 5)
  
  # Plot
  ggplot2::ggplot(df, ggplot2::aes(lib, count, col = virus)) +
    ggplot2::geom_point(position = "jitter", size = 2) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::theme_classic(15) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust=1, vjust=1),
      panel.background = ggplot2::element_rect(fill = "grey90"),
      panel.grid.major.x = ggplot2::element_line(color = "white", size = 8),
      panel.grid.major.y = ggplot2::element_line(colour = "black")
    ) +
    ggplot2::labs(x = NULL)
  fn <- paste0(out_dir, "viral-counts-sum", ".", "png")
  ggplot2::ggsave(fn, width = 15, height = 5)
  
  message("Done.")
}

if (sys.nframe() == 0) {
  main()
}
