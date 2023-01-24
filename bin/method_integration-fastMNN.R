"
Integration by fastMNN

Usage:
    method_integration-fastMNN.R [options] <file>
    
Options:
    -h --help     Show this screen.
" -> doc

# Source functions
suppressMessages({
  for (i in list.files("bin/functions", full.names = TRUE)) {
    source(i)
  }
})

#' Run fastMNN
#' 
#' @param ds Dataset of class 'SingleCellExperiment'
#' 
#' @returns Dataset of class 'SingleCellExperiment'
#' 
run_fastMNN <- function(ds = NULL, d = 30, k = 10, assay = "X", 
                        hvg="highly_variable") {
  
  # Run fastMNN
  features <- rownames(ds)[which(ds@rowRanges@elementMetadata[, hvg] == TRUE)]
  sce <- batchelor::fastMNN(
    ds, batch=ds$batch, subset.row=features, d=d, k=k, assay.type=assay
    )
  
  emb <- data.frame(sce@int_colData$reducedDims$corrected)
  names(emb) <- paste0("X_mnn", 1:ncol(emb))
  
  return(emb)
}


#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  args <- docopt::docopt(doc)
  
  in_file <- args[["<file>"]]
  out_file <- in_file
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_file, obsm = FALSE, layers = "cp10k")
  
  # Run method -----------------------------------------------------------------
  message("Running fastMNN...")
  
  assays <- c("cp10k")
  hvg_method <- ds@metadata$hvg_keys$keys
  dims <- c(10,30,50)
  neighbors <- c(10)
  for (assay in assays) {
    for (hvg in hvg_method) {
      for (n in dims) {
        for (k in neighbors) {
          print(paste(n, "dimensions of", k, "neighbors from", hvg, "features",
                      "from assay", assay))
          key <- paste0(assay, "_", hvg, "_fastMNN", "_n", n, "k", k)
          ds@int_colData$reducedDims[[key]] <- run_fastMNN(ds, assay = assay, 
                                                           d=n, k=k, hvg=hvg)
        }
      }
    }
  }
  
  # Write output ---------------------------------------------------------------
  
  # Add embedding
  write_obsm_h5ad(ds, out_file)
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}