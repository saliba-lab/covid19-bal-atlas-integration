"
Quality control

Usage:
    quality-control.R --out-file=<path> [options] <file>
    
Options:
    -h --help             Show this screen.
    -o --out-file=<path>  Output file
" -> doc

# Source functions
suppressMessages({
  source("bin/_functions.R")
})

#' Load data
#' 
#' @param input Path to AnnData file
#' 
#' @returns SingleCellExperiment
read_data <- function(input) {
  
  # Fetch components -----------------------------------------------------------
  ds <- read_layer_h5ad(input,"X")
  rn <- rownames(ds)
  cn <- colnames(ds)
  dimnames(ds) <- list(NULL, NULL)
  obs <- read_slot_h5ad(input, "obs")
  var <- read_slot_h5ad(input, "var")
  scov <- read_slot_h5ad(input, "SCoV2_counts")
  
  # Create SCE -----------------------------------------------------------------
  ds <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ds),
    colData = obs, rowData = var
  )
  
  scov <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(scov))
    )
  SingleCellExperiment::altExp(ds, "SCoV2") <- scov
  
  return(ds)
  
}

#' Add quality metrics
#' 
#' @param input SingleCellExperiment
#' 
#' @returns The function output
add_quality_metrics <- function(input) {
  
  
  
}

#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  args <- docopt::docopt(doc)
  
  file_in <- args[["<file>"]]
  file_out <- args[["--out-file"]]
  
  # Read data ------------------------------------------------------------------
  
  message("Reading data from '", file_in, "'...")
  ds <- read_data(file_in)
  message("Read data:")
  print(ds)
  
  message("Writing output to '", out_file, "'...")
  # write_output(output, out_file)
  message("Done!")
}

if (sys.nframe() == 0) {
  #main()
  print("Yay")
}