# Download viral genomes

# Set paths
file <- "docs/download/viral-assembly_genome-urls.csv"
out.fa <- "docs/genomes/viral.fasta"
out.gtf <- "docs/genomes/viral.gtf"
dir <- "docs/genomes/"
dir.create(dir, recursive = TRUE)

# Load data
list <- read.csv(file)

# Download files ---------------------------------------------------------------
message("Downloading FASTA & GTF files...")

for (i in 1:nrow(list)) {
  name <- list[i, "name"]
  type <- list[i, "type"]
  id <- list[i, "refseq"]
  key <- paste0(dir, name, ".", type)
  print(key)
  
  url <- list[i, "url"]
  
  download.file(url, key)
}

# Combined FASTA files ---------------------------------------------------------
message("Combining FASTA files...")

call <- paste0("cat ", dir, "*.fasta", " > ", out.fa)
system(call)

# Create GTF for all virus genomes ---------------------------------------------
message("Creating GTF file...")

# TODO: Create gtf from gff3 files
# So far done manually...

message("Done.")
