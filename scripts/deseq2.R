
# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

txi_file <- args[1]
samplesheet_file <- args[2]
output_file <- args[3]

.libPaths(c("/N/slate/gkarutur/Rlibs", .libPaths()))
library(DESeq2)



# Load tximport object
txi <- readRDS(txi_file)

# Read sample metadata
sample_info <- read.csv(
    samplesheet_file,
    stringsAsFactors = FALSE
)

# Keep only the columns needed for this comparison
sample_info <- sample_info[, c("sample", "condition")]

# Set sample names as row names
rownames(sample_info) <- sample_info$sample

# Define reference condition
sample_info$condition <- factor(
    sample_info$condition,
    levels = c("PBS", "ZIKV")
)

# Match metadata order to tximport count matrix
sample_info <- sample_info[
    colnames(txi$counts),
    ,
    drop = FALSE
]

# Check sample matching
if (!identical(rownames(sample_info), colnames(txi$counts))) {
    stop("Sample names do not match between tximport and metadata")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(
    txi,
    colData = sample_info,
    design = ~ condition
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

# Run differential expression analysis
dds <- DESeq(dds)

# Extract ZIKV versus PBS results
res <- results(
    dds,
    contrast = c("condition", "ZIKV", "PBS")
)

# Convert results to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Sort by adjusted p-value
res_df <- res_df[
    order(res_df$padj, na.last = TRUE),
]

# Save results
write.csv(
    res_df,
    output_file,
    row.names = FALSE
)




# Print summary
cat("DESeq2 analysis completed\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes tested:", nrow(dds), "\n")
cat(
    "Significant genes (padj < 0.05):",
    sum(res_df$padj < 0.05, na.rm = TRUE),
    "\n"
)
cat(
    "Upregulated:",
    sum(res_df$padj < 0.05 &
        res_df$log2FoldChange > 0, na.rm = TRUE),
    "\n"
)
cat(
    "Downregulated:",
    sum(res_df$padj < 0.05 &
        res_df$log2FoldChange < 0, na.rm = TRUE),
    "\n"
)


# Create volcano plot
library(ggplot2)

## Define significance categories
res_df$significance <- "Not significant"

res_df$significance[
    !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange > 0
] <- "Upregulated"

res_df$significance[
    !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange < 0
] <- "Downregulated"

# Avoid infinite values when adjusted p-value is zero
res_df$padj_plot <- pmax(
    res_df$padj,
    .Machine$double.xmin,
    na.rm = FALSE
)

# Generate plot
volcano <- ggplot(
    res_df,
    aes(x = log2FoldChange, y = -log10(padj_plot),
        color = significance)
) +
    geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    scale_color_manual(values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not significant" = "grey"
    )) +
    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed"
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +
    labs(
        title = "ZIKV vs PBS: Differential Gene Expression",
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value",
        color = "Gene category"
    ) +
    theme_classic()

# Save plot
ggsave(
    "volcano_plot.png",
    plot = volcano,
    width = 8,
    height = 6,
    dpi = 300
)

cat("Volcano plot saved: volcano_plot.png\n")

