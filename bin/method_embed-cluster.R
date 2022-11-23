"
Calculate embeddings & clusterings

Usage:
    method_embed-cluster.R [options] <file>
    
Options:
    -h --help             Show this screen.
    -o --out-dir=<path>   Output directory
" -> doc

# Source functions
suppressMessages({
  source("bin/_functions.R")
})

#' Read AnnData
#' 
#' @param input Path to AnnData file
#' 
#' @returns SingleCellExperiment
read_h5ad <- function(input = NULL) {
  
  stopifnot(!is.null(input))
  
  # Fetch components
  ds <- read_layer_h5ad(input,"X")
  dimnames(ds) <- list(NULL, NULL)
  obs <- read_slot_h5ad(input, "obs")
  var <- read_slot_h5ad(input, "var")
  
  # Create SCE
  ds <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = ds),
    colData = obs, rowData = var
  )
  
  # Add obsm
  file <- rhdf5::h5ls(input)
  names <- file$name[file$group == "/obsm"]
  names <- names[!names == "SCoV2_counts"]
  
  for (name in names) {
    print(paste("Adding", name))
    mat <- as.matrix(t(
      rhdf5::h5read(input, paste0("obsm/", name))
    ))
    rownames(mat) <- rownames(obs)
    colnames(mat) <- paste0("component_", 1:ncol(mat))
    
    SingleCellExperiment::reducedDim(ds, name) <- mat
  }
  
  return(ds)
}

#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  args <- docopt::docopt(doc)
  
  in_file <- args[["<file>"]]
  out_file <- in_file
  out_dir <- args[["--out-dir"]]
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  if (!endsWith(out_dir, "/")) {
    out_dir <- paste0(out_dir, "/")
  }
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_file)
  
  # Calculate UMAP embedding & clusters ----------------------------------------
  set.seed(42)
  dimreds <- c("X_pca", "X_mnn", "X_scVI")
  components <- c(10, 20, 30, 40, 50)
  settings <- data.frame(
    row.names = dimreds,
    min_dist = c(0.5, 0.5, 0.3)
  )
  algorithm <- "leiden"
  
  for (i in dimreds) {
    for (n in components) {
      # UMAP embedding
      emb <- paste0(i, "_", n, "_", "umap")
      ds <- scater::runUMAP(ds, name=emb, dimred=i, n_dimred=n, 
                            min_dist=settings[i, "min_dist"])
      # Clustering
      clust <- paste0(i, "_", n, "_", algorithm)
      ds@metadata[[clust]] <- cluster_cells(ds, emb=i, alg=algorithm, dims = n,
                                          res_from=0.1, res_to=4, res_by=0.1)
      # Plot
      plot_clusters_embedding(ds, clusters = clust, embedding = emb, 
                              write = TRUE, out_dir = out_dir, nrow = 4)
    }
  }
  
  # Write output ---------------------------------------------------------------
  
  # View file
  file <- rhdf5::h5ls(out_file)
  file[file$group == "/obsm", ]
  
  # Add embeddings
  for (i in names(ds@int_colData$reducedDims)) {
    message(paste("Adding", i))
    key <- paste0("obsm/", i)
    
    if (i %in% file$name) {
      rhdf5::h5delete(out_file, key)
    }
    
    obj <- as.matrix(t(ds@int_colData$reducedDims[[i]]))
    rhdf5::h5write(obj = obj, file = out_file, name = key)
  }
  
  # Add clusters
  for (i in names(ds@metadata)) {
    message(paste("Adding", i))
    key <- paste0("obsm/", i)
    
    if (i %in% file$name) {
      rhdf5::h5delete(out_file, key)
    }
    
    obj <- as.matrix(t(ds@metadata[[i]]))
    rhdf5::h5write(obj = obj, file = out_file, name = key)
  }
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}