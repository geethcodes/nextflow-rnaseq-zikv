args <- commandArgs(trailingOnly = TRUE)

tx2gene_file <- args[1]
output_file <- args[2]
quant_dirs <- args[-c(1, 2)]

library(tximport)

files <- file.path(quant_dirs, "quant.sf")

names(files) <- basename(quant_dirs)

if (!all(file.exists(files))) {
    stop("Missing Salmon quantification files")
}

tx2gene <- read.delim(
    tx2gene_file,
    header = FALSE,
    col.names = c("TXNAME", "GENEID")
)

txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "no"
)

saveRDS(txi, output_file)

cat("TXIMPORT completed\n")
cat("Samples:", length(files), "\n")
cat("Genes:", nrow(txi$counts), "\n")
