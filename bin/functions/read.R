#' Get layer from AnnData
#' 
#' @param file Path to h5ad file of type AnnData
#' @param assay Assay name (default: X)
#' 
#' @returns Sparse matrix
#' 
read_layer_h5ad <- function(file = NULL, name = "X") {
  
  stopifnot(
    endsWith(file, "h5ad"),
    name %in% rhdf5::h5ls(file)$name
  )
  
  if (name != "X") {
    name <- paste0("layers/", name)
  }
  
  X <- Matrix::sparseMatrix(
    i = rhdf5::h5read(file, paste0(name, "/indices")) + 1,
    p = rhdf5::h5read(file, paste0(name, "/indptr")),
    x = as.integer(rhdf5::h5read(file, paste0(name, "/data"))),
    dims = c(
      length(rhdf5::h5read(file, "var/_index")), 
      length(rhdf5::h5read(file, "obs/_index"))
    )
  )
  
  return(X)
}


#' Get slot from AnnData
#' 
#' @param file Path to h5ad file of type AnnData
#' @param slot Slot name (default: obs)
#' 
#' @return Data.frame
#' 
read_slot_h5ad <- function(
    file = NULL,
    slot = "obs"
) {
  
  stopifnot(
    endsWith(file, "h5ad"),
    slot %in% rhdf5::h5ls(file)$name
  )
  
  index <- match(slot, rhdf5::h5ls(file)$name)
  path <- paste0(rhdf5::h5ls(file)$group[index], "/", slot)
  data <- rhdf5::h5read(file, path, read.attributes = TRUE)
  
  # Convert lists into factors
  for (i in names(data)) {
    if (class(data[[i]]) == "list") {
      print(paste("Converting", i, "to factor"))
      index <- data[[i]]$codes + 2
      lvls <- c(NA, data[[i]]$categories)
      names(lvls) <- sort(unique(index))
      data[[i]] <- factor(lvls[index], lvls)
    }
  }
  
  # Convert to encoding type
  df <- data.frame(
    row.names = data[[attr(data, "_index")]]
  )
  for (i in attr(data, "column-order")) {
    df[[i]] <- data[[i]]
  }
  
  return(df)
}


#' Read obsm from AnnData
#' 
#' @param file Path to h5ad file of type AnnData
#' 
#' @returns List of matrices
#' 
read_obsm_h5ad <- function(file) {
  
  stopifnot(
    endsWith(file, "h5ad")
  )
  
  data <- rhdf5::h5read(file, "obsm", read.attributes = TRUE)
  
  obsm <- list()
  for (i in names(data)) {
    
    if (all(class(data[[i]]) == "list")) {
      obsm[[i]] <- data.frame(
        row.names = data[[i]][["_index"]]
      )
      for (j in names(data[[i]])) {
        obsm[[i]][[j]] <- data[[i]][[j]]
      }
      obsm[[i]][["_index"]] <- NULL
    }
    
    if (any(class(data[[i]]) == "matrix")) {
      obsm[[i]] <- t(data[[i]])
    }
  }
  attributes(obsm) <- attributes(data)
  
  return(obsm)
}


#' Read AnnData
#' 
#' @param input Path to AnnData file
#' 
#' @returns SingleCellExperiment
read_h5ad <- function(file = NULL,
                      layers = FALSE,
                      obsm = TRUE,
                      uns = TRUE
                      ) {
  
  stopifnot(
    endsWith(file, "h5ad")
    )
  
  # File overview
  fs <- rhdf5::h5ls(file)
  
  # Select main layer
  all_layers <- fs$name[fs$group == "/layers"]
  if (!layers %in% c(TRUE, FALSE, all_layers)) {
    stop(paste(layers, "not available in", file))
  }
  if (layers %in% c(TRUE, FALSE)) {
    main_layer <- "X"
  } else {
    main_layer <- layers
    layers <- FALSE
  }
  
  # Create SCE
  ds <- SingleCellExperiment::SingleCellExperiment(
    assays = list(X = read_layer_h5ad(file, main_layer)),
    colData = read_slot_h5ad(file, "obs"), 
    rowData = read_slot_h5ad(file, "var")
  )
  
  # Add layers
  if (layers) {
    assays <- fs$name[fs$group == "/layers"]
    for (i in assays) {
      ds@assays@data[[i]] <- read_layer_h5ad(file, i)
    }
  }
  
  # Add obsm
  if (obsm) {
    obsm <- read_obsm_h5ad(file)
    for (i in names(obsm)) {
      SingleCellExperiment::reducedDim(ds, i) <- obsm[[i]]
    }
  }
  
  # Add uns
  if (uns) {
    ds@metadata <- rhdf5::h5read(file, "uns")
  }
  
  return(ds)
}
