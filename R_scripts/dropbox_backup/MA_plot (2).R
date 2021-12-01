library("DESeq2")
library("geneplotter")

dds <- DESeqDataSet(se[, 1:6], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[1:6]
dds_DESeq <- DESeq(dds)
res_GV <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
plotMA(res_GV, ylim = c(-2, 2), cex = 0.65)
savepng("MA_plot_GV_KO_vs_WT", width = 1000, asp = 0.75)

dds <- DESeqDataSet(se[, 7:12], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[7:12]
dds_DESeq <- DESeq(dds)
res_MII <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
plotMA(res_MII, ylim = c(-2, 2), cex = 0.65)
savepng("MA_plot_MII_KO_vs_WT", width = 1000, asp = 0.75)

dds <- DESeqDataSet(se[, 13:18], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[13:18]
dds_DESeq <- DESeq(dds)
res_1C <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
plotMA(res_1C, ylim = c(-2, 2), cex = 0.65)
savepng("MA_plot_oneC_KO_vs_WT", width = 1000, asp = 0.75)

plotMA(res_GV, ylim = c(-2, 2), cex = 0.65)
