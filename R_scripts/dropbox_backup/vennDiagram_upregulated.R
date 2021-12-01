library("VennDiagram")
library("geneplotter")
setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/DESeq_analysis_Petr/Venn_diagrams")

genes_GV <- read.csv("GV.csv")
genes_MII <- read.csv("MII.csv")
genes_1C <- read.csv("1C.csv")

genes_GV <- genes_GV[which(!is.na(genes_GV$entrezID)), ]
genes_GV <- genes_GV[, 1:8]
upregulated_genes_GV <- genes_GV[genes_GV$log2FoldChange > 0, ]

genes_MII <- genes_MII[which(!is.na(genes_MII$entrezID)), ]
genes_MII <- genes_MII[, 1:8]
upregulated_genes_MII <- genes_MII[genes_MII$log2FoldChange > 0, ]

genes_1C <- genes_1C[which(!is.na(genes_1C$entrezID)), ]
genes_1C <- genes_1C[, 1:8]
upregulated_genes_1C <- genes_1C[genes_1C$log2FoldChange > 0, ]

# plot three samples
library("VennDiagram")
library("geneplotter")
grid.newpage()
venn.plot <- venn.diagram(list(upregulated_genes_GV$entrezID, upregulated_genes_MII$entrezID, upregulated_genes_1C$entrezID), 
                          NULL, 
                          fill = c("blue", "red", "green"), 
                          alpha = c(0.5, 0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("GV", "MII", "1C"),
                          main = "Upregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("all_upregulated", width = 1000, asp = 1)

# plot two samples
grid.newpage()
venn.plot <- venn.diagram(list(upregulated_genes_GV$entrezID, upregulated_genes_MII$entrezID), 
                          NULL, 
                          fill = c("blue", "red"), 
                          alpha = c(0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("GV", "MII"),
                          main = "Upregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("GV_MII_upregulated", width = 1000, asp = 1)

grid.newpage()
venn.plot <- venn.diagram(list(upregulated_genes_GV$entrezID, upregulated_genes_1C$entrezID), 
                          NULL, 
                          fill = c("blue", "green"), 
                          alpha = c(0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("GV", "1C"),
                          main = "Upregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("GV_1C_upregulated", width = 1000, asp = 1)

grid.newpage()
venn.plot <- venn.diagram(list(upregulated_genes_MII$entrezID, upregulated_genes_1C$entrezID), 
                          NULL, 
                          fill = c("red", "green"), 
                          alpha = c(0.5, 0.5), 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 3,
                          category.names = c("MII", "1C"),
                          main = "Upregulated genes", 
                          main.cex = 3)
grid.draw(venn.plot)
savepng("MII_1C_upregulated", width = 1000, asp = 1)

# tables 

# GV
upregulated_GV_only <- upregulated_genes_GV[-which(upregulated_genes_GV$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_MII$entrezID)), ]
upregulated_GV_only <- upregulated_GV_only[-which(upregulated_GV_only$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_1C$entrezID)), ]
write.csv(upregulated_GV_only, "upregulated_GV_only.csv")

upregulated_GV_MII_only <- upregulated_genes_GV[which(upregulated_genes_GV$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_MII$entrezID)), ]
upregulated_GV_MII_only <- upregulated_GV_MII_only[-which(upregulated_GV_MII_only$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_1C$entrezID)), ]
write.csv(upregulated_GV_MII_only, "upregulated_GV_MII_only.csv")

upregulated_GV_1C_only <- upregulated_genes_GV[which(upregulated_genes_GV$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_1C$entrezID)), ]
upregulated_GV_1C_only <- upregulated_GV_1C_only[-which(upregulated_GV_1C_only$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_MII$entrezID)), ]
write.csv(upregulated_GV_1C_only, "upregulated_GV_1C_only.csv")

upregulated_GV_all <- upregulated_genes_GV[which(upregulated_genes_GV$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_MII$entrezID)), ]
upregulated_GV_all <- upregulated_GV_all[which(upregulated_GV_all$entrezID %in% intersect(upregulated_GV_all$entrezID, upregulated_genes_1C$entrezID)), ]
write.csv(upregulated_GV_all, "upregulated_GV_all.csv")

# MII
upregulated_MII_only <- upregulated_genes_MII[-which(upregulated_genes_MII$entrezID %in% intersect(upregulated_genes_MII$entrezID, upregulated_genes_GV$entrezID)), ]
upregulated_MII_only <- upregulated_MII_only[-which(upregulated_MII_only$entrezID %in% intersect(upregulated_genes_MII$entrezID, upregulated_genes_1C$entrezID)), ]
write.csv(upregulated_MII_only, "upregulated_MII_only.csv")

upregulated_MII_GV_only <- upregulated_genes_MII[which(upregulated_genes_MII$entrezID %in% intersect(upregulated_genes_MII$entrezID, upregulated_genes_GV$entrezID)), ]
upregulated_MII_GV_only <- upregulated_MII_GV_only[-which(upregulated_MII_GV_only$entrezID %in% intersect(upregulated_genes_GV$entrezID, upregulated_genes_1C$entrezID)), ]
write.csv(upregulated_MII_GV_only, "upregulated_MII_GV_only.csv")

upregulated_MII_1C_only <- upregulated_genes_MII[which(upregulated_genes_MII$entrezID %in% intersect(upregulated_genes_MII$entrezID, upregulated_genes_1C$entrezID)), ]
upregulated_MII_1C_only <- upregulated_MII_1C_only[-which(upregulated_MII_1C_only$entrezID %in% intersect(upregulated_genes_MII$entrezID, upregulated_genes_GV$entrezID)), ]
write.csv(upregulated_MII_1C_only, "upregulated_MII_1C_only.csv")

upregulated_MII_all <- upregulated_genes_MII[which(upregulated_genes_MII$entrezID %in% intersect(upregulated_genes_MII$entrezID, upregulated_genes_GV$entrezID)), ]
upregulated_MII_all <- upregulated_MII_all[which(upregulated_MII_all$entrezID %in% intersect(upregulated_MII_all$entrezID, upregulated_genes_1C$entrezID)), ]
write.csv(upregulated_MII_all, "upregulated_MII_all.csv")


# 1C
upregulated_one_1C_only <- upregulated_genes_1C[-which(upregulated_genes_1C$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_GV$entrezID)), ]
upregulated_one_1C_only <- upregulated_one_1C_only[-which(upregulated_one_1C_only$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_MII$entrezID)), ]
write.csv(upregulated_one_1C_only, "upregulated_1C_only.csv")

upregulated_one_1C_GV_only <- upregulated_genes_1C[which(upregulated_genes_1C$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_GV$entrezID)), ]
upregulated_one_1C_GV_only <- upregulated_one_1C_GV_only[-which(upregulated_one_1C_GV_only$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_MII$entrezID)), ]
write.csv(upregulated_one_1C_GV_only, "upregulated_1C_GV_only.csv")

upregulated_one_1C_MII_only <- upregulated_genes_1C[which(upregulated_genes_1C$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_MII$entrezID)), ]
upregulated_one_1C_MII_only <- upregulated_one_1C_MII_only[-which(upregulated_one_1C_MII_only$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_GV$entrezID)), ]
write.csv(upregulated_one_1C_MII_only, "upregulated_1C_MII_only.csv")

upregulated_one_1C_all <- upregulated_genes_1C[which(upregulated_genes_1C$entrezID %in% intersect(upregulated_genes_1C$entrezID, upregulated_genes_GV$entrezID)), ]
upregulated_one_1C_all <- upregulated_one_1C_all[which(upregulated_one_1C_all$entrezID %in% intersect(upregulated_one_1C_all$entrezID, upregulated_genes_MII$entrezID)), ]
write.csv(upregulated_one_1C_all, "upregulated_1C_all.csv")
