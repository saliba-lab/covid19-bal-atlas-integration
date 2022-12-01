# Download viral genomes

# Set paths
file <- "docs/viral-assembly_genome-urls.csv"
dir <- "data/genomes/viral/"
dir.create(dir)

# Load data
list <- read.csv(file)

# Download files
for (i in 1:nrow(list)) {
  name <- list[i, "name"]
  type <- list[i, "type"]
  id <- list[i, "refseq"]
  key <- paste0(dir, name, ".", type)
  print(key)
  
  url <- list[i, "url"]
  
  download.file(url, key)
}
