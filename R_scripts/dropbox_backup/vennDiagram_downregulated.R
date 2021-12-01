library("VennDiagram")
library("geneplotter")
setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/DESeq_analysis_Petr/Venn_diagrams")

genes_GV <- read.csv("GV.csv")
genes_MII <- read.csv("MII.csv")
genes_1C <- read.csv("1C.csv")

genes_GV <- genes_GV[which(!is.na(genes_GV$entrezID)), ]
genes_GV <- genes_GV[, 1:8]
downregulated_genes_GV <- genes_GV[genes_GV$log2FoldChange < 0, ]

genes_MII <- genes_MII[which(!is.na(genes_MII$entrezID)), ]
genes_MII <- genes_MII[, 1:8]
downregulated_genes_MII <- genes_MII[genes_MII$log2FoldChange < 0, ]

genes_1C <- genes_1C[which(!is.na(genes_1C$entrezID)), ]
genes_1C <- genes_1C[, 1:8]
downregulated_genes_1C <- genes_1C[genes_1C$log2FoldChange < 0, ]


# plot three samples
grid.newpage()
venn.plot <- venn.diagram(list(downregulated_genes_GV$entrezID, downregulated_genes_MII$entrezID, downregulated_genes_1C$entrezID), 
                          NULL, 
                          fill = c("blue", "red", "green"), 
                          alpha = c(0.5, 0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("GV", "MII", "1C"),
                          main = "Downregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("all_downregulated", width = 1000, asp = 1)

# plot two samples
grid.newpage()
venn.plot <- venn.diagram(list(downregulated_genes_GV$entrezID, downregulated_genes_MII$entrezID), 
                          NULL, 
                          fill = c("blue", "red"), 
                          alpha = c(0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("GV", "MII"),
                          main = "Downregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("GV_MII_downregulated", width = 1000, asp = 1)

grid.newpage()
venn.plot <- venn.diagram(list(downregulated_genes_GV$entrezID, downregulated_genes_1C$entrezID), 
                          NULL, 
                          fill = c("blue", "green"), 
                          alpha = c(0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("GV", "1C"),
                          main = "Downregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("GV_1C_downregulated", width = 1000, asp = 1)

grid.newpage()
venn.plot <- venn.diagram(list(downregulated_genes_MII$entrezID, downregulated_genes_1C$entrezID), 
                          NULL, 
                          fill = c("red", "green"), 
                          alpha = c(0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("MII", "1C"),
                          main = "Downregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("MII_1C_downregulated", width = 1000, asp = 1)


# tables
# GV
downregulated_GV_only <- downregulated_genes_GV[-which(downregulated_genes_GV$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_MII$entrezID)), ]
downregulated_GV_only <- downregulated_GV_only[-which(downregulated_GV_only$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_1C$entrezID)), ]
write.csv(downregulated_GV_only, "downregulated_GV_only.csv")

downregulated_GV_MII_only <- downregulated_genes_GV[which(downregulated_genes_GV$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_MII$entrezID)), ]
downregulated_GV_MII_only <- downregulated_GV_MII_only[-which(downregulated_GV_MII_only$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_1C$entrezID)), ]
write.csv(downregulated_GV_MII_only, "downregulated_GV_MII_only.csv")

downregulated_GV_1C_only <- downregulated_genes_GV[which(downregulated_genes_GV$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_1C$entrezID)), ]
downregulated_GV_1C_only <- downregulated_GV_1C_only[-which(downregulated_GV_1C_only$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_MII$entrezID)), ]
write.csv(downregulated_GV_1C_only, "downregulated_GV_1C_only.csv")

downregulated_GV_all <- downregulated_genes_GV[which(downregulated_genes_GV$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_MII$entrezID)), ]
downregulated_GV_all <- downregulated_GV_all[which(downregulated_GV_all$entrezID %in% intersect(downregulated_GV_all$entrezID, downregulated_genes_1C$entrezID)), ]
write.csv(downregulated_GV_all, "downregulated_GV_all.csv")

# MII
downregulated_MII_only <- downregulated_genes_MII[-which(downregulated_genes_MII$entrezID %in% intersect(downregulated_genes_MII$entrezID, downregulated_genes_GV$entrezID)), ]
downregulated_MII_only <- downregulated_MII_only[-which(downregulated_MII_only$entrezID %in% intersect(downregulated_genes_MII$entrezID, downregulated_genes_1C$entrezID)), ]
write.csv(downregulated_MII_only, "downregulated_MII_only.csv")

downregulated_MII_GV_only <- downregulated_genes_MII[which(downregulated_genes_MII$entrezID %in% intersect(downregulated_genes_MII$entrezID, downregulated_genes_GV$entrezID)), ]
downregulated_MII_GV_only <- downregulated_MII_GV_only[-which(downregulated_MII_GV_only$entrezID %in% intersect(downregulated_genes_GV$entrezID, downregulated_genes_1C$entrezID)), ]
write.csv(downregulated_MII_GV_only, "downregulated_MII_GV_only.csv")

downregulated_MII_1C_only <- downregulated_genes_MII[which(downregulated_genes_MII$entrezID %in% intersect(downregulated_genes_MII$entrezID, downregulated_genes_1C$entrezID)), ]
downregulated_MII_1C_only <- downregulated_MII_1C_only[-which(downregulated_MII_1C_only$entrezID %in% intersect(downregulated_genes_MII$entrezID, downregulated_genes_GV$entrezID)), ]
write.csv(downregulated_MII_1C_only, "downregulated_MII_1C_only.csv")

downregulated_MII_all <- downregulated_genes_MII[which(downregulated_genes_MII$entrezID %in% intersect(downregulated_genes_MII$entrezID, downregulated_genes_GV$entrezID)), ]
downregulated_MII_all <- downregulated_MII_all[which(downregulated_MII_all$entrezID %in% intersect(downregulated_MII_all$entrezID, downregulated_genes_1C$entrezID)), ]
write.csv(downregulated_MII_all, "downregulated_MII_all.csv")

# 1C
downregulated_one_1C_only <- downregulated_genes_1C[-which(downregulated_genes_1C$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_GV$entrezID)), ]
downregulated_one_1C_only <- downregulated_one_1C_only[-which(downregulated_one_1C_only$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_MII$entrezID)), ]
write.csv(downregulated_one_1C_only, "downregulated_1C_only.csv")

downregulated_one_1C_GV_only <- downregulated_genes_1C[which(downregulated_genes_1C$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_GV$entrezID)), ]
downregulated_one_1C_GV_only <- downregulated_one_1C_GV_only[-which(downregulated_one_1C_GV_only$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_MII$entrezID)), ]
write.csv(downregulated_one_1C_GV_only, "downregulated_1C_GV_only.csv")

downregulated_one_1C_MII_only <- downregulated_genes_1C[which(downregulated_genes_1C$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_MII$entrezID)), ]
downregulated_one_1C_MII_only <- downregulated_one_1C_MII_only[-which(downregulated_one_1C_MII_only$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_GV$entrezID)), ]
write.csv(downregulated_one_1C_MII_only, "downregulated_1C_MII_only.csv")

downregulated_one_1C_all <- downregulated_genes_1C[which(downregulated_genes_1C$entrezID %in% intersect(downregulated_genes_1C$entrezID, downregulated_genes_GV$entrezID)), ]
downregulated_one_1C_all <- downregulated_one_1C_all[which(downregulated_one_1C_all$entrezID %in% intersect(downregulated_one_1C_all$entrezID, downregulated_genes_MII$entrezID)), ]
write.csv(downregulated_one_1C_all, "downregulated_1C_all.csv")
