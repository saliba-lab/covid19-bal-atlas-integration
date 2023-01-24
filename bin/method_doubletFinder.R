"
Doublet detection

Usage:
    analysis_annotation-celltype.R [options] <file>
    
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
  out_file <- stringr::str_replace(in_file, ".h5ad", "_doublets.tsv")
  
  # Read data ------------------------------------------------------------------
  message("Reading data...")
  
  ds <- read_h5ad(in_file, layers = "counts")
  gc()
  
  # Find doublets --------------------------------------------------------------
  message("Running scDblFinder...")
  
  names(ds@assays) <- "counts"
  ds <- scDblFinder::scDblFinder(ds, samples="libname", nfeatures = 2000)
  
  index <- names(ds@colData)
  data <- ds@colData[, index[stringr::str_detect(index, "scDbl")]]
  
  # Write data -----------------------------------------------------------------
  message(paste("Writing data to", out_file))
  
  write.table(data, out_file)
  
  message("Done!")
  }

if (sys.nframe() == 0) {
  main()
}