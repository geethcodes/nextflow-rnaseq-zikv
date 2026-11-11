suppressPackageStartupMessages({ library(DESeq2); library(org.Mm.eg.db) })
dds <- readRDS("results/deseq2/dds.rds")
res <- results(dds, contrast = c("condition","ZIKV","PBS"))
sig_all <- rownames(res)[which(res$padj < 0.05)]

# 1. Assign fetal sex from Y-chromosome genes
ygenes <- mapIds(org.Mm.eg.db, c("Ddx3y","Eif2s3y","Kdm5d","Uty"), "ENSEMBL", "SYMBOL")
ycnt <- colSums(counts(dds, normalized = TRUE)[intersect(ygenes, rownames(dds)), , drop = FALSE])
dds$sex <- factor(ifelse(ycnt > 0.2 * max(ycnt), "M", "F"))
print(data.frame(sample = colnames(dds), condition = dds$condition, Y_counts = round(ycnt), sex = dds$sex))
print(table(dds$condition, dds$sex))

# 2. Sex as covariate: design ~ sex + condition
d_sex <- dds; design(d_sex) <- ~ sex + condition
d_sex <- DESeq(d_sex, quiet = TRUE)
r_sex <- results(d_sex, contrast = c("condition","ZIKV","PBS"))
sig_sex <- rownames(r_sex)[which(r_sex$padj < 0.05)]

# 3. Without the outlier PBS_2
d2 <- dds[, colnames(dds) != "PBS_2"]; d2$condition <- droplevels(d2$condition)
d2 <- DESeq(d2, quiet = TRUE)
r2 <- results(d2, contrast = c("condition","ZIKV","PBS"))
sig2 <- rownames(r2)[which(r2$padj < 0.05)]

lfc_r <- function(a, b) round(cor(a$log2FoldChange, b$log2FoldChange[match(rownames(a), rownames(b))], use = "complete.obs"), 3)
cat("\n                         sig   overlap_with_original   LFC_r\n")
cat(sprintf("Original ~condition    %5d\n", length(sig_all)))
cat(sprintf("~sex + condition       %5d   %5d                   %s\n", length(sig_sex), length(intersect(sig_all, sig_sex)), lfc_r(res, r_sex)))
cat(sprintf("Without PBS_2          %5d   %5d                   %s\n", length(sig2), length(intersect(sig_all, sig2)), lfc_r(res, r2)))

# 4. Technical signature check: fold change vs transcript length
len <- rowMeans(assays(dds)[["avgTxLength"]])
cat("\nSpearman LFC vs transcript length:", round(cor(res$log2FoldChange, len, method = "spearman", use = "complete.obs"), 3), "\n")

write.csv(as.data.frame(r_sex), "results/deseq2/DESeq2_results_sex_adjusted.csv")
write.csv(as.data.frame(r2),    "results/deseq2/DESeq2_results_without_PBS2.csv")
