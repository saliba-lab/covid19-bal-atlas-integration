"
Quality control for the combined dataset

Usage:
    analysis-combined-qc.R [options]
    
Options:
    -h --help     Show this screen.
" -> doc

# Global variables -------------------------------------------------------------

# Input params
file <- "docs/overview.xlsx"
url <- "https://nubes.helmholtz-berlin.de/s/LnJn5z8wB2o2NrJ"
url <- paste0(url, "/download")

# Output params
fastq_path <- "docs/fastq-path.csv"
sheet_path <- "docs/samplesheets/"
plot_dir <- paste0("analysis/combined/overview/")
dir.create(plot_dir, recursive = TRUE)

# Download data ----------------------------------------------------------------

# Download
download.file(url, file)

# Prepare fastq paths ----------------------------------------------------------

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

# Subset to columns included in samplesheet
index <- match(c("path_1", "run", "path_2", "flowcell", "sample"), names(data))
fastqs <- data[, index]

# Combine into single path
collapse <- function (x) {stringr::str_c(x, collapse = "/")}
fastqs <- apply(fastqs, 1, collapse)

# Remove whitespaces
fastqs <- stringr::str_remove(fastqs, "\\s")

# Save file
writeLines(fastqs, fastq_path)

# Prepare run samplesheets -----------------------------------------------------

# Subset to columns included in samplesheet
index <- match(c("lane", "sample", "index"), names(data))
sheets <- data[, index]
names(sheets) <- stringr::str_to_title(names(sheets))

# Split samplesheet by run
sheets <- split(sheets, data$run)

# Save samplesheets as CSV
dir.create(sheet_path, recursive = TRUE)
for (i in names(sheets)) {
  csv <- paste0(sheet_path, i, ".csv")
  write.csv(sheets[[i]], csv, row.names = FALSE)
}

# Plot sample/patient overview -------------------------------------------------

# Create output directory
dir.create(dest, recursive = TRUE)

# Load patient data
patients <- readxl::read_excel(file, "patient")

# Add patient data to samples
index <- match(data$patient, patients$id)
data$age <- patients$age[index]
data$sex <- patients$sex[index]

# Create sample id
data$type_sample <- paste0(data$type, "-", data$sample, " (", data$patient, ")")

# Re-order data
data <- data[order(data$dpso_value), ]
data$type_sample <- factor(data$type_sample, unique(data$type_sample))

sind <- c("Male" = "\u2642", "Female" = "\u2640")

# Plot sample data
ggplot2::ggplot(data, ggplot2::aes(dpso_value, type_sample, fill = age)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::geom_point(
    ggplot2::aes(x = -5, color = sex), 
    shape = sind[data$sex], size = 5
  ) +
  ggplot2::theme_classic(20) +
  ggplot2::labs(x = "Time post symptom onset", y = "Patient") +
  ggplot2::guides(
    fill = ggplot2::guide_legend(override.aes = list(size = 0))
  ) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(),
    panel.grid.minor.x = ggplot2::element_line(linetype = "dashed")
  )

fn <- paste0(plot_dir, "patients.png")
ggplot2::ggsave(fn, width = 12, height = 6)
