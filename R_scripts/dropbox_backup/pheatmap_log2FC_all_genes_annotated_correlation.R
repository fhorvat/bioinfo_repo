library("gplots")
library("pheatmap")
library("NMF")
library("geneplotter")
library("dynamicTreeCut")
library("WGCNA")
library("RColorBrewer")

results_all_no_NA <- results_all[rowSums(is.na(results_all)) < 1, ]

dist_mat <- as.dist(1 - abs(cor(t(results_all_no_NA))))
hc.rows <- hclust(dist_mat, method = "average")

# Dynamic Tree Cut
maxCoreScatter <- 0.93
minGap <- (1 - maxCoreScatter) * 3/4
dynamicCut <- cutreeDynamic(hc.rows,
                            method = "hybrid", 
                            distM = as.matrix(dist_mat), 
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
savepng("dtc_abs_correlation_average", width = 1000, asp = 0.75)

ann_row <- data.frame(cluster = factor(dynamicCut),
                      row.names = rownames(results_all_no_NA))
ann_colors <- list(cluster = c("1" = "turquoise", "2" = "blue", "3" = "brown", 
                               "4" = "yellow", "5" = "green", "6" = "red", 
                               "7" = "black", "8" = "pink", "9" = "magenta", 
                               "10" = "purple", "11" = "greenyellow"))

bk <- unique(c(seq(-2, -0.2, length = 10), seq(-0.2, 0.2, length = 20), seq(0.2, 2, length = 10)))
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)

pheatmap_output <- pheatmap(results_all_no_NA, 
                            clustering_method = "average",
                            # clustering_distance_rows = dist_mat,
                            col = hmcols, 
                            breaks = bk,
                            fontsize = 15,
                            row_clusters = T,
                            treeheight_row = 300, 
                            show_rownames = F,
                            cluster_cols = F,
                            annotation_row = ann_row,
                            annotation_colors = ann_colors)
savepng("pheatmap_abs_correlation_average", width = 1000, asp = 0.75)