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
  
  # Extract flowcell name from run and remove leading character (A & B)
  data$flowcell <- stringr::str_split(data$run, "_", simplify = TRUE)[, 4]
  data$flowcell <- gsub("^\\w", "", data$flowcell)
  
  # Add path elements
  data$path_1 <- "data/raw/fastq"
  data$path_2 <- "outs/fastq_path"
  
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
  
  # Manage exception for run 221014
  # Problem was the use of single-index (i5) for HTOs with
  # dual indices (SI-TT-) for gene expression
  run <- "221014_A00643_0567_AHYMN3DSX3"
  sheet <- sheets[[run]]
  
  key <- paste0(run, "_HTO")
  sheets[[key]] <- sheet[stringr::str_which(sheet$Index,"SI-TT",negate=TRUE), ]
  sheets[[run]] <- sheet[stringr::str_which(sheet$Index, "SI-TT"), ]
  
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
                   "run", "path_2", "flowcell", "index"), names(data))
  libs <- data[, index]
  libs$wd <- getwd()
  libs <- split(libs, libs$libname)
  for (i in names(libs)) {
    lib <- libs[[i]]
    csv <- data.frame(
      fastqs = stringr::str_c(
        lib[, c("wd", "path_1", "run", "path_2", "flowcell", "libname")], 
        collapse = "/"),
      sample = i,
      library_type = "Gene Expression"
    )
    if (stringr::str_detect(lib$index, "SI-TT")) {
      csv$fastqs <- stringr::str_remove(csv$fastqs, paste0("/", csv$sample))
    }
    if (!is.na(lib$`TotalSeq-A`)) {
      lib$run <- paste0(lib$run, "_HTO")
      csv <- rbind(csv, c(
        stringr::str_c(
          lib[, c("path_1", "run", "path_2")], 
          collapse = "/"),
        paste0(i, "-HTO"), "Antibody Capture"
      ))
    }
    # Add csv to list
    libs[[i]] <- csv
  }
  
  return(libs)
}

#' Create sample overview
#' 
#' @param data Data.frame of sequenced libraries
#' 
create_sample_overview <- function(data) {
  
  # Combine library, sample & patient data
  patients <- readxl::read_excel(file, "patients")
  samples  <- readxl::read_excel(file, "samples")
  
  # Subset to sequenced samples
  sind <- stringr::str_split(data$sample, ",")
  names(sind) <- 1:length(sind)
  for (i in names(sind)) {
    if (length(sind[[i]]) > 1) {
      sind[[i]] <- stringr::str_split(sind[[i]], ":", simplify = TRUE)[,2]
    }
  }
  sind <- unlist(sind)
  samples <- samples[samples$sample %in% sind, ]
  samples <- samples[samples$cohort %in% c("A", "B", "C", "Recovery"), ]
  
  # Calculate DPSO
  samples$onset <- patients[["Hospital-admission"]][match(
    samples$patient, patients$id
  )]
  ind <- is.na(samples$dpso)
  dpso <- samples$date - samples$onset
  samples$dpso[ind] <- dpso[ind]
  
  # Calculate duration to endpoint
  samples$endpoint_date <- patients$endpoint_date[match(
    samples$patient, patients$id
  )]
  samples$endpoint_dpso <- samples$endpoint_date - samples$onset
  samples$endpoint_dpso[is.na(samples$endpoint_dpso)] <- -10
  
  # Add patient information
  ind <- match(samples$patient, patients$id)
  samples$sex <- patients$sex[ind]
  samples$endpoint <- patients$endpoint[ind]
  
  # Order
  ind <- dplyr::summarise(
    dplyr::group_by(samples, patient), dpso = max(dpso), cohort = cohort,
    endpoint = endpoint
  )
  ind <- unique(ind)
  ind <- ind$patient[order(ind$cohort, ind$endpoint, ind$dpso)]
  samples$patient <- factor(as.character(samples$patient), ind)
  
  # Plot
  plot <- ggplot2::ggplot(samples, ggplot2::aes(dpso, patient)) +
    ggplot2::geom_point(ggplot2::aes(fill = cohort), size = 3, shape = 21) +
    ggplot2::scale_fill_manual(values = color$cohort) +
    ggplot2::geom_point(
      ggplot2::aes(x = -5, color = sex), 
      shape = shape$sex[samples$sex], size = 4
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = endpoint_dpso, shape = endpoint), size = 4
    ) +
    ggplot2::scale_shape_manual(values = shape$endpoint) +
    ggplot2::theme_classic(15) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      panel.background = ggplot2::element_rect(fill = "grey90"),
      panel.grid.major.y = ggplot2::element_line(color = "white", size = 5),
      panel.grid.major.x = ggplot2::element_line(size = 1, color = "grey50"),
      panel.grid.minor.x = ggplot2::element_line(size = .5, color = "grey50")
    ) +
    ggplot2::labs(x = "Days post symptom onset", y = "Patient")
  
  return(plot)
}

#' The main script function
main <- function() {
  
  # Variables ------------------------------------------------------------------
  file <- "docs/overview.xlsx"
  url <- "https://nubes.helmholtz-berlin.de/s/LnJn5z8wB2o2NrJ"
  url <- paste0(url, "/download")
  
  mkfastq_samplesheets <- "data/samplesheets/mkfastq/"
  count_library_csv <- "data/samplesheets/count/"
  hto_indices <- "docs/HTO-indices.csv"
  plot_dir <- paste0("analysis/BCB/overview/")
  
  # Create directories
  dir.create(plot_dir, recursive = TRUE)
  dir.create(mkfastq_samplesheets, recursive = TRUE)
  dir.create(count_library_csv, recursive = TRUE)
  
  # Colors & shapes ------------------------------------------------------------
  
  # Color
  color <- list()
  color$cohort <- c(
    "A" = "indianred",
    "B" = "orange",
    "C" = "orangered",
    "Recovery" = "purple"
  )
  
  # Shape
  shape <- list()
  shape$sex <- c(
    "Male" = "\u2642", 
    "Female" = "\u2640"
    )
  shape$endpoint <- c(
    death = "\u271D",
    release = "\u00BB",
    unknown = "\u003F"
  )
  
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
  
  # Sample overview ------------------------------------------------------------
  
  create_sample_overview(data)
  
  # Save
  fn <- paste0(plot_dir, "sample-overview.png")
  ggplot2::ggsave(fn, width = 12, height = 6)
  
  message("Done!")
}

if (sys.nframe() == 0) {
  main()
}

