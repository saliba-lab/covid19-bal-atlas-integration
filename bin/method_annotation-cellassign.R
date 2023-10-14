"
Cell type annotation

Usage:
    analysis_annotation-cellassign.R [options] <file>
    
Options:
    -h --help             Show this screen.
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
  
  out_file_ann <- stringr::str_replace(out_file, ".h5ad", "_cellassign.csv")
  out_dir <- stringr::str_replace(
    stringr::str_replace(in_file, "data", "analysis"), ".h5ad", "/"
  )
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_file, layers = "counts")
  
  # Marker genes
  xls_url <- "https://nubes.helmholtz-berlin.de/s/id9rrKNeditKMmF"
  xls_url <- paste0(xls_url, "/", "download")
  xls <- "docs/celltype/markers.xlsx"
  download.file(xls_url, xls)
  markers <- list()
  for (i in readxl::excel_sheets(xls)) {
    markers[[i]] <- readxl::read_excel(xls, sheet = i)
  }
  
  # Run CellAssign -------------------------------------------------------------
  message("Running CellAssign...")
  
  mrk <- markers$level_3
  mrk <- split(mrk$gene, mrk$type)
  mrk <- cellassign::marker_list_to_mat(mrk)
  
  genes <- rownames(mrk)
  genes <- genes[sparseMatrixStats::rowSums2(ds[genes, ]@assays@data$X)>0]
  
  cells <- colnames(ds)
  cells <- cells[sparseMatrixStats::colSums2(ds[genes, ]@assays@data$X)>0]
  
  mrk <- mrk[genes, ]
  label <- cellassign::cellassign(ds[genes, cells], mrk, sce_assay = "X")
  
  cellassign::cellassign(
    cellassign::example_sce, cellassign::example_rho
  )
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}
