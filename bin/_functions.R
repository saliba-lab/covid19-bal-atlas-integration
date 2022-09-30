#' Get assay from anndata h5ad file
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
read_assay_h5ad <- function(
    file = NULL,
    assay = "X"
    ) {
  
  stopifnot(
    !is.null(file),
    assay %in% rhdf5::h5ls(file)$name
  )
  
  row_ind <- rhdf5::h5read(file, "var/_index")
  col_ind <- rhdf5::h5read(file, "obs/_index")
  X <- Matrix::sparseMatrix(
    i = rhdf5::h5read(file, "X/indices") + 1,
    p = rhdf5::h5read(file, "X/indptr"),
    x = as.numeric(rhdf5::h5read(file, "X/data"))
  )
  
  rownames(X) <- row_ind
  colnames(X) <- col_ind
  
  return(X)
  
}

#' Get 'var' from anndata h5ad file
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
  
  data <- rhdf5::h5read(file, slot)
  attrs <- rhdf5::h5readAttributes(file, slot)
  
  # Convert lists into factors
  for (i in names(data)) {
    if (class(data[[i]]) == "list") {
      print(paste("Converting", i, "to factor"))
      index <- data[[i]]$codes + 1
      lvls <- data[[i]]$categories
      data[[i]] <- factor(lvls[index], lvls)
    }
  }
  
  # TODO: Check conversion to factors in obs
  
  # Convert to encoding type
  df <- data.frame(
    row.names = data[[attrs[["_index"]]]]
  )
  for (i in attrs[["column-order"]]) {
    df[[i]] <- data[[i]]
  }
  
  return(df)
  
}
