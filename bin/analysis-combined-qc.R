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

# Global variables -------------------------------------------------------------
file_input_adata <- "data/combined.h5ad"
file_input_aggr_csv <- "docs/combined_aggr.csv"
file_input_overview <- "docs/overview.xlsx"

plot_output <- "analysis/combined/qc/"
dir.create(plot_output, recursive = TRUE)

file_output_obs <- "docs/combined_coldata.csv"

# Load data --------------------------------------------------------------------

ds <- read_layer_h5ad(file_input_adata)
obs <- read_slot_h5ad(file_input_adata, "obs")
scov <- read_slot_h5ad(file_input_adata, "SCoV2_counts")

# Convert barcode tag to sample name
data <- read.csv(file_input_aggr_csv)
index <- data$library_id
names(index) <- rownames(data)

obs$tag <- stringr::str_split(colnames(ds), "-", simplify = TRUE)[, 2]
obs$sample <- index[obs$tag]

# Add sample metadata ----------------------------------------------------------

readxl::excel_sheets(file_input_overview)
data <- readxl::read_excel(file_input_overview, "sample")

names(data)
cols_keep <- c("type", "patient", "dpso_value", "status", "cohort")

obs <- cbind(obs, data[match(obs$sample, data$sample), cols_keep])

# Add cell quality metrics -----------------------------------------------------

# Library size
obs$libsize <- Matrix::colSums(ds)

# Features
obs$nfeatures <- Matrix::colSums(ds > 0)

# Percent mitochrondrial counts
index <- which(stringr::str_detect(rownames(ds), "^MT-"))
obs$percent.mt <- round(Matrix::colSums(ds[index, ]) / obs$libsize, 3) * 100

# Percent viral count
obs$scov2_counts <- rowSums(scov)
obs$percent.scov2 <- round(obs$scov2_counts / obs$libsize, 3) * 100

# Assign cell quality ----------------------------------------------------------

# Set thresholds
lin.t <- data.frame(
  row.names  = c("min", "max"),
  libsize    = c(0, Inf),
  features   = c(200, Inf),
  percent.mt = c(0, 25)
)
loess.t <- -0.15

for (t_type in c("linear", "adaptive", "loess")) {
  
  print(t_type)
  
  # Pull data
  data <- obs[, c("libsize", "nfeatures", "percent.mt")]
  names(data) <- c("libsize", "features", "percent.mt")
  
  if (t_type == "linear") {
    # Make assignment
    data$quality <- dplyr::between(data$libsize, lin.t$libsize[1], lin.t$libsize[2]) &
      dplyr::between(data$features, lin.t$features[1], lin.t$features[2]) & 
      dplyr::between(data$percent.mt, lin.t$percent.mt[1], lin.t$percent.mt[2])
  } 
  
  if (t_type == "adaptive") {
    df <- scuttle::perCellQCFilters(
      data, sum.field = "libsize", detected.field = "features", 
      sub.fields = "percent.mt"
    )
    data$quality <- df$discard == FALSE
  }
  
  if (t_type == "loess") {
    
    # Remove quality assignment
    data$quality <- NULL
    
    # Loess fit
    fit <- loess(features ~ libsize, log10(data), span = 1)
    fit <- as.data.frame(fit[c("x", "y", "fitted", "residuals")])
    
    # Model fit
    ggplot2::ggplot(fit, ggplot2::aes(libsize, y)) +
      ggplot2::geom_point(shape = 1, alpha = .5, size = .5) +
      ggplot2::geom_line(ggplot2::aes(y = fitted), col = "blue")
    fn <- paste0(plot_output, t_type, "_", "model-fit", ".", "png")
    ggplot2::ggsave(fn, width = 6, height = 6)
    
    
    # Loess threshold
    ggplot2::ggplot(fit, ggplot2::aes(libsize, residuals)) +
      ggplot2::geom_point(shape = 1, alpha = .5, size = .5) +
      ggplot2::geom_line(ggplot2::aes(y = 0), col = "blue") +
      ggplot2::geom_line(ggplot2::aes(y = loess.t), col = "red")
    fn <- paste0(plot_output, t_type, "_", "model-thresh", ".", "png")
    ggplot2::ggsave(fn, width = 6, height = 6)
    
    # Store quality assignments
    data$quality <- fit$residuals > loess.t
  }
  
  data$quality <- factor(
    c("TRUE" = "good", "FALSE" = "bad")[as.character(data$quality)]
  )
  
  # Violin plots
  df <- tidyr::gather(data, "metric", "value", -quality)
  ggplot2::ggplot(
    df, ggplot2::aes(metric, value, col = quality)
  ) +
    ggplot2::geom_point(position = "jitter", shape = 16, size = .1) +
    ggplot2::geom_violin(
      scale = "width", col = "black", size = 1, fill = NA, 
      draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2
    ) +
    ggplot2::facet_wrap(~metric, scales = "free") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::theme_classic(20) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    )
  
  fn <- paste0(plot_output, t_type, "_", "qc-violin", ".", "png")
  ggplot2::ggsave(fn, width = 12, height = 6)
  
  # Correlation scatterplots
  df <- tidyr::gather(data, "metric", "value", -libsize, -quality)
  ggplot2::ggplot(
    df, ggplot2::aes(libsize, value, col = quality)
  ) +
    ggplot2::geom_point(shape = 16, size = .1) +
    ggplot2::facet_wrap(~metric, scales = "free") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::theme_classic(20) +
    ggplot2::labs(x = "Library size", y = NULL) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    )
  
  fn <- paste0(plot_output, t_type, "_", "qc-scatter", ".", "png")
  ggplot2::ggsave(fn, width = 12, height = 6)
  
  # Add quality to ds
  key <- paste0("qc.", t_type)
  obs[[key]] <- data$quality
  
  # Compare across samples -----------------------------------------------------
  
  # Add patient IDs
  data$patient <- obs$patient
  
  # Libsize - features
  ggplot2::ggplot(data, ggplot2::aes(libsize, features)) +
    ggplot2::geom_point(ggplot2::aes(col = quality), size = .1) +
    ggplot2::facet_wrap(~patient) +
    ggplot2::theme_classic(10) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    )
  
  fn <- paste0(
    plot_output, t_type, "_", "patient", "_", "libsize-features", ".", "png"
    )
  ggplot2::ggsave(fn, width = 8, height = 8)
  
  ggplot2::ggplot(data, ggplot2::aes(libsize, percent.mt)) +
    ggplot2::geom_point(ggplot2::aes(col = quality), size = .1) +
    ggplot2::facet_wrap(~patient) +
    ggplot2::theme_classic(10) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    )
  
  fn <- paste0(
    plot_output, t_type, "_", "patient", "_", "libsize-percMT", ".", "png"
    )
  ggplot2::ggsave(fn, width = 8, height = 8)
  
}

# Write output files -----------------------------------------------------------

# Add random selection
bcs <- split(rownames(obs), obs$sample)
bcs <- unlist(lapply(bcs, function(x) {sample(x, 600)}))
obs$random <- rownames(obs) %in% bcs

# Metadata
write.csv(obs, file_output_obs)
