library("dynamicTreeCut")
library("DESeq2")
library("WGCNA")
library("pheatmap")
library("geneplotter")
library("RColorBrewer")

dds <- DESeqDataSet(se[, 1:6], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[1:6]
dds_DESeq <- DESeq(dds)
res_GV <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))

dds <- DESeqDataSet(se[, 7:12], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[7:12]
dds_DESeq <- DESeq(dds)
res_MII <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))

dds <- DESeqDataSet(se[, 13:18], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[13:18]
dds_DESeq <- DESeq(dds)
res_1C <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))

significant_rownames <- unique(c(rownames(subset(res_GV, padj < 0.1)), 
                                 rownames(subset(res_MII, padj < 0.1)),
                                 rownames(subset(res_1C, padj < 0.1))))
res_GV_significant <- res_GV[which(rownames(res_GV) %in% significant_rownames), ]
res_MII_significant <- res_MII[which(rownames(res_MII) %in% significant_rownames), ]
res_1C_significant <- res_1C[which(rownames(res_1C) %in% significant_rownames), ]
res_all_significant <- cbind(res_GV_significant$"log2FoldChange", 
                             res_MII_significant$"log2FoldChange",
                             res_1C_significant$"log2FoldChange")
colnames(res_all_significant) <- c("GV", "MII", "1C")
rownames(res_all_significant) <- rownames(res_GV_significant)

# pheatmap
bk <- unique(c(seq(-2, -0.2, length = 10), seq(-0.2, 0.2, length = 30), seq(0.2, 2, length = 10)))
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)
# hmcols <- colorRampPalette(c("red", "white", "blue"))(length(bk) - 1)

pheatmap(res_all_significant, 
         clustering_method = "complete",
         col = hmcols, 
         breaks = bk,
         cluster_row = T, 
         cluster_cols = F,
         treeheight_row = 300, 
         show_rownames = F,
         main = "log2FC DESeq2 significant complete method")
savepng("pheatmap_log2FC_significant_genes_complete", width = 1000)

res_all_significant_dist <- as.matrix(res_all_significant)
res_all_significant_dist <- dist(res_all_significant_dist)
hc.rows <- hclust(res_all_significant_dist, method = "ward.D")
maxCoreScatter <- 0.90
minGap <- (1 - maxCoreScatter) * 3/4
dynamicCut <- cutreeDynamic(hc.rows,
                            method = "hybrid", 
                            distM = as.matrix(res_all_significant_dist), 
                            deepSplit = 2, 
                            maxCoreScatter = maxCoreScatter, 
                            minGap = minGap, 
                            maxAbsCoreScatter = NULL, 
                            minAbsGap = NULL)
cut2colour <- labels2colors(dynamicCut) 
plotDendroAndColors(hc.rows, 
                    cut2colour, 
                    "Dynamic Tree Cut", 
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05) 
savepng("significant_genes_all_values_dtc_90", width = 1000, asp = 0.75)
write.csv(data.frame(res_all_significant, dynamicCut), "significant_genes_all_values_dtc_90")

