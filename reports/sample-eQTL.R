"
Overview of the combined dataset

Usage:
    overview_sample.R [options]
    
Options:
    -h --help     Show this screen.
" -> doc

# Source functions
suppressMessages({
  source("bin/_functions.R")
})

#' Plot sex by XIST expression
#' 
#' @param ds Dataset of type 'SingleCellExperiment'
#' 
#' @returns plot
#' 
plot_sex <- function(ds, assay="logcounts") {
  
  data <- ds@colData[order(ds$dpso_value, ds$patient), ]
  data <- data.frame(
    patient = ds$id,
    count = ds@assays@data[[assay]][grep("XIST", rownames(ds)), ]
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(count, patient)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_classic(20) +
    ggplot2::labs(title = "XIST") +
    ggplot2::scale_x_continuous(trans = "log10")
  
  return(plot)
}

#' The main script function
main <- function() {
  
  # Read data
  file <- "data/combined-core.h5ad"
  ds <- read_h5ad(file)
  ds@assays@data$counts <- read_layer_h5ad(file, "counts")
  
  # Add sample ID
  ds$id <- stringr::str_c("BAL-", ds$sample, " (", ds$patient, ")")
  ds$id <- factor(ds$id, unique(ds$id[order(ds$dpso_value)]))
  
  # Plot
  plot_sex(ds, assay = "counts")
  fn <- "analysis/combined/overview//patient-XIST.png"
  ggplot2::ggsave(fn, height = 6, width = 5)
  
}

if (sys.nframe() == 0) {
  main()
}