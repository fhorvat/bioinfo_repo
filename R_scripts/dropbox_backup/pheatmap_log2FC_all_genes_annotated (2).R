library("pheatmap")
library("geneplotter")
library("dynamicTreeCut")
library("WGCNA")
library("RColorBrewer")

results_all_no_NA <- results_all[rowSums(is.na(results_all)) < 1, ]
results_all_output_distance <- as.matrix(results_all_no_NA)
results_all_output_distance <- dist(results_all_output_distance)
hc.rows <- hclust(results_all_output_distance, method = "ward.D")

# Dynamic Tree Cut
maxCoreScatter <- 0.93
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
savepng("dtc_clusters_all_genes_log2FC_wardD", width = 1000)
all_genes_log2FC_clusters <- data.frame(results_all_no_NA, cluster = dynamicCut)
write.csv(all_genes_log2FC_clusters, "all_genes_log2FC_clusters_wardD.csv")

ann_row <- data.frame(cluster = factor(dynamicCut), row.names = rownames(results_all_no_NA))
cluster_colors <- labels2colors(1:40) 
names(cluster_colors) <- 1:40 
cluster_colors <- cluster_colors[1:length(unique(cut2colour))]
ann_colors <- list(cluster = cluster_colors)

bk <- unique(c(seq(-2, -0.2, length = 10), seq(-0.2, 0.2, length = 30), seq(0.2, 2, length = 10)))
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)
pheatmap(results_all_no_NA, 
         clustering_method = "ward.D",
         col = hmcols, 
         breaks = bk,
         fontsize = 15,
         row_clusters = T,
         treeheight_row = 300, 
         show_rownames = F,
         annotation_row = ann_row, 
         annotation_colors = ann_colors,
         cluster_cols = F)
savepng("heatmap_clusters_all_genes_log2FC_wardD", width = 1000)

# library("gplots")
# heatmap.2_output <- heatmap.2(as.matrix(results_all_clusters[, 1:3]),
#                               hclustfun = function(d) hclust(d, method = "ward"),
#                               Rowv = T, 
#                               Colv = F,
#                               breaks = bk, 
#                               symbreaks = T, 
#                               trace = "none", 
#                               RowSideColors = cut2colour, 
#                               col = mcols,
#                               labRow = "")
# savepng("heatmap2_output_ryb", width = 1000)

# library("NMF")
# aheatmap_output <- aheatmap(results_all_clusters[, 1:3], 
#                             hclustfun = "ward", 
#                             color = hmcols, 
#                             breaks = bk, 
#                             treeheight = 300, 
#                             Rowv = T, 
#                             Colv = NA, 
#                             labRow = NULL,
#                             annColors = ann_colors,
#                             annRow = ann_row)
# savepng("aheatmap_output_ryb", width = 1000)

