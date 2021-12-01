library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("geneplotter")

rpkm_new <- read.csv("rpkm_new.csv", row.names = 1)

rpkm_new[rpkm_new == 0] <- NA
rpkm_new <- do.call(data.frame, lapply(rpkm_new, function(x) replace(x, is.infinite(x), NA)))
rownames(rpkm_new) <- rownames(rpkm_new)
rpkm_new <- rpkm_new[rowSums(is.na(rpkm_new)) == 0, ]
rpkm_new <- log2(rpkm_new)
colnames(rpkm_new) <- c(colnames(rpkm_new)[1:12], paste0("oneC_KO", 1:3), paste0("oneC_WT", 1:3))

bk <- seq(range(rpkm_new)[1], range(rpkm_new)[2], length = 50)
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)

# dist_mat_rows <- as.dist(1 - cor(t(rpkm_new)))
# dist_mat_cols <- as.dist(1 - cor((rpkm_new)))

pheatmap(rpkm_new,
         clustering_method = "ward.D2",
         # clustering_distance_rows = dist_mat_rows,
         # clustering_distance_cols = dist_mat_cols,
         main = "log2 RPKM all samples ward.D2 method",
         col = hmcols, 
         breaks = bk,
         cluster_row = T, 
         cluster_cols = T,
         treeheight_row = 100, 
         treeheight_col	= 100,
         show_rownames = F)
savepng("pheatmap_log2_RPKM_all_samples_wardD2", width = 1000)


