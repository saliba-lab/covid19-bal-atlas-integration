"
Download BAL dataset from Wendisch et al. 2021 (PMID:  34914922)

Usage:
    dataset_wdm21A.R [options]
    
Options:
    -h --help     Show this screen.
" -> doc

# Source functions
suppressMessages({
  source("bin/_functions.R")
})

#' Get COVID-19 BAL macrophages
#' 
#' @returns Dataset of class 'SingleCellExperiment'
get_macrophages <- function() {
  
  # Variables
  url <- "https://nubes.helmholtz-berlin.de/s/cF5DXLZzN5EkAYZ"
  url <- paste0(url, "/download")
  tmp_file <- tempfile()
  
  # Read data
  download.file(url, tmp_file)
  ds <- readRDS(tmp_file)
  file.remove(tmp_file)
  
  # Convert to SingleCellExperiment
  ds <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = ds@assays$RNA@counts,
      logcounts = ds@assays$RNA@data
      ),
    colData = ds@meta.data,
    rowData = ds@assays$RNA@meta.features,
    reducedDims = list(
      "pca" = ds@reductions$pca@cell.embeddings,
      "umap" = ds@reductions$umap@cell.embeddings
    ),
    altExps = list(
      "SCoV2" = SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = ds@assays$SCoV2@counts)
      )
    )
  )
  
  return(ds)
}

#' The main script function
main <- function() {
  
  # Variables
  out_file <- "data/wdm21A_macrophages.h5ad"
  
  # Load data
  ds <- get_macrophages()
  
  # Write data
  write_h5ad(ds, out_file)
  
}

if (sys.nframe() == 0) {
  main()
}