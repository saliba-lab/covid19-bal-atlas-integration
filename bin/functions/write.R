#' Write matrix into obsm
#' 
#' @param file Anndata file (h5ad)
#' 
write_obsm_h5ad <- function(ds = NULL, file = NULL) {
  
  stopifnot(
    endsWith(file, "h5ad")
  )
  
  fs <- rhdf5::h5ls(file)
  file_obsm <- fs$name[fs$group == "/obsm"]
  
  message("Adding embeddings...")
  for (i in names(ds@int_colData$reducedDims)) {
    if (i %in% file_obsm) {
      message(paste("Embedding", i, "already present. Skipping..."))
    } else {
      message(paste("Adding", i, "to obsm"))
      emb <- t(ds@int_colData$reducedDims[[i]])
      rhdf5::h5write(emb, file, paste0("obsm/", i))
    }
  }
}

#' Add columns to obs
#' 
#' @param obs Data.frame with observations
#' @param file Anndata file (h5ad)
#' 
add_obs_h5ad <- function(data = NULL, file = NULL) {
  
  stopifnot(
    class(data) == "data.frame",
    endsWith(file, "h5ad")
  )
  
  # Look into file
  rhdf5::h5ls(file)
  fs <- rhdf5::h5ls(file)
  attrs <- rhdf5::h5readAttributes(file, "obs/sample_id")
  
  # Add sample
  rhdf5::h5delete(file, "obs/sample")
  v <- factor(data$sample)
  vl <- list(
    "categories" = levels(v),
    "codes" = array(match(v, levels(v)) - 1)
  )
  attr(vl, "encoding-type") <- "categorical"
  attr(vl, "encoding-version") <- "0.2.2"
  rhdf5::h5write(vl, file, "obs/sample", write.attributes = TRUE)
  
  # Edit column-order
  f <- rhdf5::H5Fopen(file)
  f_obs <- rhdf5::H5Gopen(f, "obs")
  name <- "column-order"
  f_obs_attr <- rhdf5::H5Aopen(f_obs, name)
  keys <- c(rhdf5::H5Aread(f_obs_attr), "sample")
  rhdf5::H5Awrite(f_obs_attr, keys)
  rhdf5::H5Aclose(f_obs_attr)
  
  rhdf5::h5closeAll()
  
  message("Done.")
}