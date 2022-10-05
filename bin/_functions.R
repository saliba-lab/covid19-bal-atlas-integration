#' Get X from anndata h5ad file
#' 
#' @param file Path to h5ad file in anndata format
#' @param assay Assay name (default: X)
#' 
#' @return Sparse matrix
#' 
#' author Oliver Dietrich
#' 
#' @export
#' 
read_layer_h5ad <- function(
    file = NULL,
    name = "X"
    ) {
  
  stopifnot(
    !is.null(file),
    name %in% rhdf5::h5ls(file)$name
  )
  
  if (name != "X") {
    name <- paste0("layers/", name)
  }
  
  row_ind <- rhdf5::h5read(file, "var/_index")
  col_ind <- rhdf5::h5read(file, "obs/_index")
  X <- Matrix::sparseMatrix(
    i = rhdf5::h5read(file, paste0(name, "/indices")) + 1,
    p = rhdf5::h5read(file, paste0(name, "/indptr")),
    x = as.numeric(rhdf5::h5read(file, paste0(name, "/data")))
  )
  
  rownames(X) <- row_ind
  colnames(X) <- col_ind
  
  return(X)
  
}


#' Get slot from anndata h5ad file
#' 
#' @param file Path to h5ad file in anndata format
#' @param slot Slot name (default: obs)
#' 
#' @return Data frame
#' 
#' author Oliver Dietrich
#' 
#' @export
#' 
read_slot_h5ad <- function(
    file = NULL,
    slot = "obs"
) {
  
  stopifnot(
    !is.null(file),
    slot %in% rhdf5::h5ls(file)$name
  )
  
  index <- match(slot, rhdf5::h5ls(file)$name)
  path <- paste0(h5ls(file)$group[index], "/", slot)
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

