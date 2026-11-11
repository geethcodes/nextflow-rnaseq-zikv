suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(ggplot2); library(pheatmap)
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
})
out  <- "results/figures"; dir.create(out,  recursive = TRUE, showWarnings = FALSE)
gdir <- "results/go";      dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS("results/deseq2/dds.rds")
vsd <- vst(dds, blind = TRUE)

# 1. PCA + sample distances
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) + geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) + theme_bw()
ggsave(file.path(out, "pca.png"), p, width = 6, height = 5, dpi = 150)
pheatmap(as.matrix(dist(t(assay(vsd)))), filename = file.path(out, "sample_distance.png"))

# 2. Annotated results
res <- results(dds, contrast = c("condition", "ZIKV", "PBS"), alpha = 0.05)
shr <- lfcShrink(dds, coef = "condition_ZIKV_vs_PBS", type = "apeglm")
df <- data.frame(gene_id = rownames(res),
                 symbol  = mapIds(org.Mm.eg.db, rownames(res), "SYMBOL", "ENSEMBL", multiVals = "first"),
                 baseMean = res$baseMean, log2FC = res$log2FoldChange,
                 log2FC_shrunk = shr$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
df <- df[order(df$padj), ]
write.csv(df, "results/deseq2/DESeq2_results_annotated.csv", row.names = FALSE)
sig <- !is.na(df$padj) & df$padj < 0.05

# 3. Volcano
df$status <- ifelse(sig & df$log2FC > 0, "Up in ZIKV", ifelse(sig & df$log2FC < 0, "Down in ZIKV", "NS"))
top <- head(df[sig, ], 15)
p <- ggplot(df, aes(log2FC_shrunk, -log10(pvalue), colour = status)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_text(data = top, aes(label = symbol), size = 2.8, vjust = -0.6, colour = "black") +
  scale_colour_manual(values = c("Up in ZIKV" = "firebrick", "Down in ZIKV" = "steelblue", NS = "grey75")) +
  labs(x = "log2 fold change (ZIKV vs PBS, apeglm)", y = "-log10 p-value", colour = NULL) + theme_bw()
ggsave(file.path(out, "volcano.png"), p, width = 7, height = 5.5, dpi = 150)

# 4. Top 50 heatmap
t50 <- head(df$gene_id[sig], 50)
m <- assay(vsd)[t50, ]
lab <- df$symbol[match(t50, df$gene_id)]; rownames(m) <- ifelse(is.na(lab), t50, lab)
pheatmap(m, scale = "row", annotation_col = as.data.frame(colData(vsd)[, "condition", drop = FALSE]),
         filename = file.path(out, "top50_heatmap.png"), fontsize_row = 6)

# 5. GO BP, up and down separately, universe = tested genes
universe <- df$gene_id[!is.na(df$padj)]
go_top <- list()
for (d in c("up", "down")) {
  g <- df$gene_id[sig & (if (d == "up") df$log2FC > 0 else df$log2FC < 0)]
  ego <- enrichGO(g, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                  ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  tab <- as.data.frame(ego)
  write.csv(tab, file.path(gdir, paste0("GO_BP_", d, ".csv")), row.names = FALSE)
  go_top[[d]] <- head(tab$Description, 10)
  if (nrow(tab) > 0)
    ggsave(file.path(out, paste0("GO_BP_", d, ".png")),
           dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP:", d, "in ZIKV")),
           width = 8, height = 6, dpi = 150)
}

# 6. One summary file to read and paste
s <- c(
  sprintf("Genes tested (non-NA padj): %d", sum(!is.na(df$padj))),
  sprintf("Significant padj<0.05: %d (up %d, down %d)", sum(sig), sum(sig & df$log2FC > 0), sum(sig & df$log2FC < 0)),
  sprintf("Significant with |shrunk LFC|>1: %d", sum(sig & abs(df$log2FC_shrunk) > 1)),
  sprintf("PC1 %d%%, PC2 %d%%", pv[1], pv[2]),
  "", "Top 20 genes:", capture.output(print(head(df[sig, c("symbol","log2FC_shrunk","padj")], 20), row.names = FALSE)),
  "", "Top GO BP up:", paste(" -", go_top$up), "", "Top GO BP down:", paste(" -", go_top$down))
writeLines(s, "results/summary.txt")
writeLines(capture.output(sessionInfo()), "results/sessionInfo_interpretation.txt")
cat(s, sep = "\n")
