
# ==========================================
# DESeq2: ZIKV versus PBS
# ==========================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop("Usage: deseq2.R txi.rds samplesheet.csv output.csv")
}

txi_file <- args[1]
samplesheet_file <- args[2]
output_file <- args[3]

library(DESeq2)
library(ggplot2)

# ------------------------------------------
# 1. Load tximport object
# ------------------------------------------

txi <- readRDS(txi_file)

# ------------------------------------------
# 2. Read sample metadata
# ------------------------------------------

sample_info <- read.csv(
    samplesheet_file,
    stringsAsFactors = FALSE
)

required_columns <- c("sample", "condition")

if (!all(required_columns %in% colnames(sample_info))) {
    stop("Samplesheet must contain sample and condition columns")
}

sample_info <- sample_info[, required_columns]

if (anyDuplicated(sample_info$sample)) {
    stop("Duplicate sample names in samplesheet")
}

if (!all(sample_info$condition %in% c("PBS", "ZIKV"))) {
    stop("Unexpected condition. Expected PBS or ZIKV.")
}

rownames(sample_info) <- sample_info$sample

sample_info$condition <- factor(
    sample_info$condition,
    levels = c("PBS", "ZIKV")
)

# ------------------------------------------
# 3. Match samples
# ------------------------------------------

count_samples <- colnames(txi$counts)
metadata_samples <- rownames(sample_info)

missing_metadata <- setdiff(count_samples, metadata_samples)
extra_metadata <- setdiff(metadata_samples, count_samples)

if (length(missing_metadata) > 0 ||
    length(extra_metadata) > 0) {

    stop(
        "Sample mismatch!\n",
        "Missing metadata: ",
        paste(missing_metadata, collapse = ", "),
        "\nExtra metadata: ",
        paste(extra_metadata, collapse = ", ")
    )
}

sample_info <- sample_info[
    count_samples,
    ,
    drop = FALSE
]

stopifnot(
    identical(rownames(sample_info), colnames(txi$counts))
)

if (any(table(sample_info$condition) < 2)) {
    stop("At least two biological replicates per condition are required")
}

# ------------------------------------------
# 4. Create DESeq2 dataset
# ------------------------------------------

dds <- DESeqDataSetFromTximport(
    txi,
    colData = sample_info,
    design = ~ condition
)

# ------------------------------------------
# 5. Filter low-count genes
# ------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 4

dds <- dds[keep, ]

cat("Genes retained after filtering:", nrow(dds), "\n")

# ------------------------------------------
# 6. Differential expression analysis
# ------------------------------------------

dds <- DESeq(dds)
saveRDS(dds, "dds.rds")

res <- results(
    dds,
    contrast = c("condition", "ZIKV", "PBS")
)

# ------------------------------------------
# 7. Export results
# ------------------------------------------

res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

res_df <- res_df[
    order(res_df$padj, na.last = TRUE),
]

write.csv(
    res_df,
    output_file,
    row.names = FALSE
)

# ------------------------------------------
# 8. Print analysis summary
# ------------------------------------------

significant <- !is.na(res_df$padj) &
    res_df$padj < 0.05

upregulated <- significant &
    res_df$log2FoldChange > 0

downregulated <- significant &
    res_df$log2FoldChange < 0

cat("\nDESeq2 analysis completed\n")
cat("Samples:", ncol(dds), "\n")
cat("Genes tested:", nrow(dds), "\n")
cat("Significant genes:", sum(significant), "\n")
cat("Upregulated:", sum(upregulated), "\n")
cat("Downregulated:", sum(downregulated), "\n")

# ------------------------------------------
# 9. Volcano plot
# ------------------------------------------

res_df$significance <- "Not significant"

res_df$significance[upregulated] <- "Upregulated"
res_df$significance[downregulated] <- "Downregulated"

res_df$padj_plot <- pmax(
    res_df$padj,
    .Machine$double.xmin
)

volcano <- ggplot(
    res_df,
    aes(
        x = log2FoldChange,
        y = -log10(padj_plot),
        color = significance
    )
) +
    geom_point(
        alpha = 0.65,
        size = 1.2,
        na.rm = TRUE
    ) +
    scale_color_manual(
        values = c(
            "Upregulated" = "red",
            "Downregulated" = "blue",
            "Not significant" = "grey"
        )
    ) +
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

ggsave(
    "volcano_plot.png",
    plot = volcano,
    width = 8,
    height = 6,
    dpi = 300
)

cat("Volcano plot saved successfully\n")

