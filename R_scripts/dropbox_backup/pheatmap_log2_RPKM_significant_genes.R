library("pheatmap")
library("geneplotter")
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

rpkm_new <- read.csv("rpkm_new.csv", row.names = 1)
rpkm_new_significant <- rpkm_new[which(rownames(rpkm_new) %in% rownames(res_all_significant)), ]

rpkm_new_significant[rpkm_new_significant == 0] <- NA
# rpkm_new_significant <- do.call(data.frame, lapply(rpkm_new_significant, function(x) replace(x, is.infinite(x), NA)))
# rownames(rpkm_new_significant) <- rownames(rpkm_new_significant)
rpkm_new_significant <- rpkm_new_significant[rowSums(is.na(rpkm_new_significant)) == 0, ]
rpkm_new_significant <- log2(rpkm_new_significant)
colnames(rpkm_new_significant) <- c(colnames(rpkm_new_significant)[1:12], 
                                    paste0("oneC_KO", 1:3), 
                                    paste0("oneC_WT", 1:3))

bk <- seq(range(rpkm_new_significant)[1], range(rpkm_new_significant)[2], length = 50)
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)

pheatmap(rpkm_new_significant,
         clustering_method = "complete",
         # clustering_distance_rows = dist_mat_rows,
         # clustering_distance_cols = dist_mat_cols,
         main = "log2 RPKM significant genes complete method",
         col = hmcols, 
         breaks = bk,
         cluster_row = T, 
         cluster_cols = T,
         treeheight_row = 100, 
         treeheight_col	= 100,
         show_rownames = F)
savepng("pheatmap_log2_RPKM_significant_genes_complete", width = 1000)

