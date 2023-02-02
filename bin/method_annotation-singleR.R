"
Cell type annotation

Usage:
    method_annotation-singleR.R [options] <ref> <query>
    
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
  
  in_ref <- args[["<ref>"]]
  in_query <- args[["<query>"]]
  
  out_file <- "..."
  
  # Read data ------------------------------------------------------------------
  message("Reading data ...")
  ds <- read_h5ad(in_query, layers = "cp10k")
  
  ref <- read_h5ad(in_ref, layers = "cp10k")
  
  # Run SingleR ----------------------------------------------------------------
  message("Running SingleR...")
  
  label <- SingleR::SingleR(query, ref,
                            assay.type.test = "X", assay.type.ref = "X"
                            )
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}
