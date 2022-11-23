"
Overview of the combined dataset

Usage:
    overview_sample.R [options]
    
Options:
    -h --help     Show this screen.
" -> doc

#' Read sample overview
#' 
#' @param file Excel file with dataset overview
#' 
#' @returns Data.frame with sample overview
#' 
read_overview_samples <- function(file) {
  
  # Read into memory
  data <- readxl::read_excel(file, "sample")
  
  # Remove empty rows
  index <- which(data$run != "")
  data <- data[index, ]
  
  # Extract flowcell name from run
  data$flowcell <- stringr::str_split(data$run, "_", simplify = TRUE)[, 4]
  
  # Remove leading character (A & B)
  data$flowcell <- gsub("^\\w", "", data$flowcell)
  
  # Add path elements
  data$path_1 <- "data/raw/fastq"
  data$path_2 <- "outs/fastq_path"
  
  return(data)
}

#' Add patient overview
#' 
#' @param file Excel file with dataset overview
#' 
#' @returns Data.frame with sample overview
#' 
add_patient_overview <- function(data, file) {
  
  # Load patient data
  patients <- readxl::read_excel(file, "patient")
  
  # Add patient data
  index <- match(data$patient, patients$id)
  data$age <- patients$age[index]
  data$sex <- patients$sex[index]
  
  # Create sample id
  data$type_sample <- paste0(data$type, "-", data$sample, " (", data$patient, ")")
  
  return(data)
}

#' Extract fastq paths
#' 
#' @param file Excel file with sample overview
#' 
#' @returns Character of FASTQ paths
#' 
extract_fastqs <- function(data) {
  
  # Subset to columns included in samplesheet
  index <- match(c("path_1", "run", "path_2", "flowcell", "sample"), names(data))
  fastqs <- data[, index]
  
  # Combine into single path
  collapse <- function (x) {stringr::str_c(x, collapse = "/")}
  fastqs <- apply(fastqs, 1, collapse)
  
  # Remove whitespaces
  fastqs <- stringr::str_remove(fastqs, "\\s")
  
  return(fastqs)
}

#' Extract sheets
#' 
#' @param data Data.frame with sample overview
#' 
#' @returns List with samplesheets
#' 
extract_sheets <- function(data) {
  
  # Subset to columns included in samplesheet
  index <- match(c("lane", "sample", "index"), names(data))
  sheets <- data[, index]
  names(sheets) <- stringr::str_to_title(names(sheets))
  
  # Split samplesheet by run
  sheets <- split(sheets, data$run)
  
  return(sheets)
}

#' Plot patient overview
#' 
#' @param data Data.frame with overview
#' 
#' @returns plot
#' 
plot_overview_patients <- function(data) {
  
  # Re-order data
  data <- data[order(data$cohort, data$dpso, data$patient), ]
  data$type_sample <- factor(data$type_sample, unique(data$type_sample))
  
  sind <- c("Male" = "\u2642", "Female" = "\u2640")
  
  # Plot sample data
  plot <- ggplot2::ggplot(data, ggplot2::aes(dpso, type_sample, fill = age)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_point(
      ggplot2::aes(x = -5, color = sex), 
      shape = sind[data$sex], size = 5
    ) +
    ggplot2::theme_classic(20) +
    ggplot2::labs(x = "Days post symptom onset", y = "Patient") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(size = 0))
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(),
      panel.grid.minor.x = ggplot2::element_line(linetype = "dashed"),
      strip.text.y = ggplot2::element_text(angle = 0)
    ) +
    ggplot2::facet_grid(cohort ~ ., scales = "free_y", space = "free_y")
  
  plot
  
  fn <- "analysis/BCB/overview/samples.png"
  ggplot2::ggsave(fn, width = 12, height = 10)
  
  return(plot)
}

#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  file <- "docs/overview.xlsx"
  url <- "https://nubes.helmholtz-berlin.de/s/LnJn5z8wB2o2NrJ"
  url <- paste0(url, "/download")
  
  fastq_path <- "docs/fastq-path.csv"
  sheet_path <- "docs/samplesheets/"
  plot_dir <- paste0("analysis/BCB/overview/")
  dir.create(plot_dir, recursive = TRUE)
  
  # Read data ------------------------------------------------------------------
  download.file(url, file)
  data <- read_overview_samples(file)
  data <- add_patient_overview(data, file)
  
  # Extract components ---------------------------------------------------------
  
  # FASTQs
  fastqs <- extract_fastqs(data)
  writeLines(fastqs, fastq_path)
  
  # Samplesheets
  sheets <- extract_sheets(data)
  dir.create(sheet_path, recursive = TRUE)
  for (i in names(sheets)) {
    csv <- paste0(sheet_path, i, ".csv")
    write.csv(sheets[[i]], csv, row.names = FALSE)
  }
  
  # Plot patient overview
  plot_overview_patients(data)
  fn <- paste0(plot_dir, "samples.png")
  ggplot2::ggsave(fn, width = 12, height = 6)
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}