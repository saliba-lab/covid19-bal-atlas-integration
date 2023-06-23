# Profibrotic macrophages in COVID-19

# Libraries ------------------------------------------------------------------

suppressMessages({
  for (i in list.files("bin/functions", full.names = TRUE)) {
    source(i)
  }
})
source("reports/sample-overview.R")

# Variables ------------------------------------------------------------------
args <- docopt::docopt(doc)

in_file <- "data/BCB/full.h5ad"
in_umap <- "analysis/BCB/full/umap-counts_seurat2000_scVI_n30l1h128.csv"
in_overview <- "docs/overview.xlsx"
in_ch1 <- list(
  "coldata" = "data/cohort-A_coldata.csv",
  "url" = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867421013830-mmc2.xlsx",
  "markers" = "docs/celltype/cohort1_macro.xlsx"
)

out_dir <- "analysis/BCB/full/cohort_A/"
if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = TRUE)}

# Colors -----------------------------------------------------------------------

color <- list()

color$celltype <- c(
  "Neutrophils"      = "darkorange2",
  "Macrophages"      = "indianred",
  "T cells"          = "deepskyblue",
  "NK cells"         = "purple",
  "Plasma cells"     = "lightseagreen",
  "Erythrocytes"     = "pink2",
  "Epithelial cells" = "goldenrod",
  "B and DC"         = "steelblue",
  "New"              = "lightgrey"
)

color$cluster <- c(
  "Monocytes"     = "#009E73",
  "Mono/Mφ"       = "#E69F00",
  "CD163/LGMN-Mφ" = "#D55E00",
  "AMφ-1"         = "#A020F0",
  "AMφ-2"         = "#0072B2",
  "Prolif. AMφ"   = "#56B4E9",
  "Other"         = "lightgrey"
)

# Read data ------------------------------------------------------------------
message("Reading data ...")
ds <- read_h5ad(in_file, layers = "cp10k")

# UMAP
ds@int_colData$reducedDims$umap <- read.csv(in_umap, row.names = "X")
key <- "umap"

# Cohort 1
ch1_cd <- read.csv(in_ch1$coldata, row.names = "X")
download.file(in_ch1$url, in_ch1$markers)
readxl::excel_sheets(in_ch1$markers)
ch1_mrk <- readxl::read_excel(in_ch1$markers, "2", skip = 2)

# Transfer labels --------------------------------------------------------------

# Subset coldata
ch1_cd <- ch1_cd[rownames(ch1_cd) %in% colnames(ds), ]

# Create index
ind <- match(rownames(ch1_cd), colnames(ds))

# Transfer cell type labels
ds$Celltype <- "New"
ds$Celltype[ind] <- ch1_cd$Celltype
ds$Celltype <- factor(ds$Celltype, levels = names(color$celltype))

# Transfer macrophage labels
ds$Macrophages <- "Other"
ds$Macrophages[ind] <- ch1_cd$Macrophages
ds$Macrophages <- factor(ds$Macrophages, levels = names(color$cluster))

# Plot
plot_embedding(ds, embedding = key, color = "Celltype", color.pal = "manual", 
               color.colors = color$celltype, color.decreasing = TRUE)
fn <- paste0(out_dir, "celltype", ".", "png")
ggplot2::ggsave(fn, width = 8, height = 6, bg = "white")

# Plot
plot_embedding(ds, embedding = key, color = "Macrophages", color.pal = "manual", 
               color.colors = color$cluster, color.decreasing = TRUE)
fn <- paste0(out_dir, "macrophages", ".", "png")
ggplot2::ggsave(fn, width = 8.5, height = 6, bg = "white")

# Plot
genes <- c("VCAN", "CCR2", "CD163", "LGMN", "MERTK", "CMKLR1", "INHBA", "FABP4")
plot_genes_embedding(ds, genes, embedding = key, nrow = 2)
fn <- paste0(out_dir, "gene-highlight", ".", "png")
ggplot2::ggsave(fn, width = 12, height = 6, bg = "white")

