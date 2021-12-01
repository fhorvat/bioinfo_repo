library("pheatmap")
library("geneplotter")
library("dynamicTreeCut")
library("WGCNA")
library("RColorBrewer")

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

results_all_output_distance <- as.matrix(res_all_significant)
results_all_output_distance <- dist(results_all_output_distance)
hc.rows <- hclust(results_all_output_distance, method = "ward.D")

# Dynamic Tree Cut
maxCoreScatter <- 0.95
minGap <- (1 - maxCoreScatter) * 3/4
dynamicCut <- cutreeDynamic(hc.rows,
                            method = "hybrid", 
                            distM = as.matrix(results_all_output_distance), 
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
savepng("dtc_clusters_significant_genes_log2FC_wardD", width = 1000)
significant_genes_log2FC_clusters <- data.frame(res_all_significant, cluster = dynamicCut)
write.csv(significant_genes_log2FC_clusters, "significant_genes_log2FC_clusters_wardD.csv")

# pheatmap with cluster annotation
ann_row <- data.frame(cluster = factor(dynamicCut), row.names = rownames(res_all_significant))
cluster_colors <- labels2colors(1:40) 
names(cluster_colors) <- 1:40 
cluster_colors <- cluster_colors[1:length(unique(cut2colour))]
ann_colors <- list(cluster = cluster_colors)

bk <- unique(c(seq(-2, -0.2, length = 10), seq(-0.2, 0.2, length = 20), seq(0.2, 2, length = 10)))
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)

pheatmap_output <- pheatmap(res_all_significant, 
                            clustering_method = "ward.D",
                            col = hmcols, 
                            breaks = bk,
                            fontsize = 15,
                            annotation_row = ann_row, 
                            annotation_colors = ann_colors,
                            row_clusters = T,
                            treeheight_row = 300, 
                            show_rownames = F,
                            cluster_cols = F)
savepng("heatmap_clusters_significant_genes_log2FC_wardD", width = 1000)
