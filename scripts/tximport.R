
# ==========================================
# TXIMPORT: Salmon transcript to gene counts
# ==========================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: tximport.R tx2gene.tsv output.rds quant_dir1 quant_dir2 ...")
}

tx2gene_file <- args[1]
output_file <- args[2]
quant_dirs <- args[-c(1, 2)]

library(tximport)

# Locate Salmon quantification files
files <- file.path(quant_dirs, "quant.sf")

# Use quantification directory names as sample identifiers
sample_names <- basename(quant_dirs)

if (anyDuplicated(sample_names)) {
    stop("Duplicate sample names detected in quantification directories")
}

names(files) <- sample_names

# Verify files exist
if (!all(file.exists(files))) {
    stop("Missing Salmon quantification files: ",
         paste(files[!file.exists(files)], collapse = ", "))
}

# Read transcript-to-gene mapping
tx2gene <- read.delim(
    tx2gene_file,
    header = FALSE,
    col.names = c("TXNAME", "GENEID"),
    stringsAsFactors = FALSE
)

# Import Salmon quantifications
txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "no"
)

# Save imported object
saveRDS(txi, output_file)

cat("TXIMPORT completed successfully\n")
cat("Samples:", length(files), "\n")
cat("Genes:", nrow(txi$counts), "\n")
cat("Sample names:", paste(names(files), collapse = ", "), "\n")

