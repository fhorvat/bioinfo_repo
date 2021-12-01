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

pheatmap_method <- function(res_all_significant, method){
  pheatmap(res_all_significant,
           clustering_method = method,
           main = paste("log2FC DESeq2 significant genes", method, "method"),
           col = hmcols, 
           breaks = bk,
           cluster_row = T, 
           cluster_cols = F,
           treeheight_row = 100, 
           treeheight_col	= 100,
           show_rownames = F)
  savepng(paste0("pheatmap_log2FC_significant_genes_", method), width = 1000)  
}
lapply(c("ward.D", "ward.D2", "complete", "average"), pheatmap_method, res_all_significant = res_all_significant)

