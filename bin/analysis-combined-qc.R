"
Quality control for the combined dataset

Usage:
    analysis-combined-qc.R [options]
    
Options:
    -h --help     Show this screen.
" -> doc

# Source functions
suppressMessages({
  source("bin/_functions.R")
})

# Global parameters ------------------------------------------------------------
file <- "data/combined.h5ad"
plot_dir <- "analysis/combined/qc/"
dir.create(plot_dir, recursive = TRUE)

# Split matrix by genome -------------------------------------------------------

# Read matrix
ds <- read_assay_h5ad(file)

# Get gene_ids
var <- read_slot_h5ad(file, "var")
index <- var$gene_ids
hs_ind <- which(stringr::str_detect(index, "ENSG"))
scov_ind <- which(stringr::str_detect(index, "ENSSASG"))

# Split matrix
scov_mat <- ds[scov_ind, ]
ds <- ds[hs_ind, ]

# Continue with QC...