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
read_libaries <- function(file) {
  
  # Read into memory
  data <- readxl::read_excel(file, "libraries")
  
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
  patients <- readxl::read_excel(file, "patients")
  
  # Add patient data
  index <- match(data$patient, patients$id)
  data$age <- patients$age[index]
  data$sex <- patients$sex[index]
  
  # Create sample id
  data$type_sample <- paste0(data$type, "-", data$sample, " (", data$patient, ")")
  
  return(data)
}

#' Create cellranger mkfastq samplesheets
#' 
#' @param data Data.frame with sample overview
#' 
#' @returns List with samplesheets
#' 
create_mkfastq_samplesheets <- function(data, hto_ind) {
  
  # Subset to columns included in samplesheet
  index <- match(c("lane", "libname", "index", "run"), names(data))
  sheets <- data[, index]
  names(sheets) <- c("Lane", "Sample", "Index", "Run")
  
  # Add TotalSeq-A HTO sequences
  for (i in sheets$Sample) {
    ind <- which(data$libname == i)
    index <- data[["TotalSeq-A"]][ind]
    if (!is.na(index)) {
      barcode <- hto_ind$barcode[hto_ind$index == index]
      row <- c(data$lane[ind], paste0(i, "-HTO"), barcode, data$run[ind])
      sheets <- rbind(sheets, row)
    }
  }
  
  # Split samplesheet by run
  sheets <- split(sheets[, -4], sheets$Run)
  
  return(sheets)
}

#' Create cellranger count library CSVs
#' 
#' @param file Excel file with sample overview
#' 
#' @returns Character of FASTQ paths
#' 
create_count_samplesheets <- function(data) {
  
  # Subset to columns included in samplesheet
  index <- match(c("libname", "TotalSeq-A", "path_1", 
                   "run", "path_2", "flowcell"), names(data))
  libs <- data[, index]
  libs <- split(libs, libs$libname)
  for (i in names(libs)) {
    lib <- libs[[i]]
    csv <- data.frame(
      fastqs = stringr::str_c(
        fastqs[[i]][, c("path_1", "run", "path_2", "flowcell", "libname")], 
        collapse = "/"),
      sample = i,
      library_type = "Gene Expression"
    )
    if (!is.na(lib$`TotalSeq-A`)) {
      csv <- rbind(csv, c(
        stringr::str_c(
          fastqs[[i]][, c("path_1", "run", "path_2", "flowcell")], 
          collapse = "/"),
        paste0(i, "-HTO"), "Antibody Capture"
      ))
    }
    # Add csv to list
    libs[[i]] <- csv
  }
  
  return(libs)
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
  file <- "docs/BCB_overview.xlsx"
  url <- "https://nubes.helmholtz-berlin.de/s/LnJn5z8wB2o2NrJ"
  url <- paste0(url, "/download")
  
  mkfastq_samplesheets <- "docs/cellranger/mkfastq/"
  count_library_csv <- "docs/cellranger/count/"
  hto_indices <- "docs/HTO-indices.csv"
  plot_dir <- paste0("analysis/BCB/overview/")
  
  # Create directories
  dir.create(plot_dir, recursive = TRUE)
  dir.create(mkfastq_samplesheets, recursive = TRUE)
  dir.create(count_library_csv, recursive = TRUE)
  
  # Read data ------------------------------------------------------------------
  download.file(url, file)
  hto_ind <- read.csv(hto_indices)
  data <- read_libaries(file)
  #data <- add_patient_overview(data, file)
  
  # Mkfastq samplesheets -------------------------------------------------------
  
  sheets <- create_mkfastq_samplesheets(data, hto_ind)
  for (i in names(sheets)) {
    csv <- paste0(mkfastq_samplesheets, i, ".csv")
    write.csv(sheets[[i]], csv, row.names = FALSE)
  }
  
  # Count samplesheets ---------------------------------------------------------
  
  sheets <- create_count_samplesheets(data)
  for (i in names(sheets)) {
    csv <- paste0(count_library_csv, i, ".csv")
    write.csv(sheets[[i]], csv, row.names = FALSE)
  }
  
  # Patient overview -----------------------------------------------------------
  
  # Detect samples across libraries
  samples <- stringr::str_split(data$sample, ",")
  names(samples) <- 1:length(samples)
  for (i in names(samples)) {
    if (length(samples[[i]]) > 1) {
      samples[[i]] <- stringr::str_split(samples[[i]], ":", simplify = TRUE)[,2]
    }
  }
  samples <- unlist(samples)
  
  # Add sample data
  pnts <- readxl::read_excel(file, "patients")
  smpl  <- readxl::read_excel(file, "samples")
  smpl$sequenced <- smpl$sample %in% unique(samples)
  
  # Plot patient overview
  plot_overview_patients(pnts)
  fn <- paste0(plot_dir, "samples.png")
  ggplot2::ggsave(fn, width = 12, height = 6)
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}