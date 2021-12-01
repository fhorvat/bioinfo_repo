library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("geneplotter")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("gage")
library("pathview")
library("PerformanceAnalytics")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/")

filenames <- file.path("/common/WORK/kristian/Projekti/Petr/Cnot6L/Mapping/bbmap/mm9", paste0("11919X", 1:18, "_sorted.bam"))
sampleTable <- read.csv("CNOT6L_sample_list_11919R_2015-10-29.csv", header = T)
sampleTable <- sampleTable[1:18, ]
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, by = "gene")

# counting
register(MulticoreParam())
se <- summarizeOverlaps(features = ebg, 
                        reads = bamfiles, 
                        mode = "Union", 
                        singleEnd = FALSE, 
                        ignore.strand = TRUE)
colData(se) <- DataFrame(sampleTable)

dds <- DESeqDataSet(se[, 7:12], design = ~Treatment.Control)
colnames(dds) <- sampleTable$ID[7:12]
dds_DESeq <- DESeq(dds)
res_MII <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_MII$sig <- "BLCK"
res_MII$sig[(!is.na(res_MII$padj) & res_MII$padj < 0.1 & res_MII$log2FoldChange > 0)] <- "RD"
res_MII$sig[(!is.na(res_MII$padj) & res_MII$padj < 0.1 & res_MII$log2FoldChange < 0)] <- "BLU"
counts_MII <- assay(se[, 7:12])
colnames(counts_MII) <- sampleTable$ID[7:12]
counts_MII <- counts_MII[res_MII$sig == "RD" | res_MII$sig == "BLU", ]
chart.Correlation(log(counts_MII + 1))

dds <- DESeqDataSet(se[, 1:6], design = ~Treatment.Control)
colnames(dds) <- sampleTable$ID[1:6]
dds_DESeq <- DESeq(dds)
res_GV <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_GV$sig <- res_MII$sig
counts_GV <- assay(se[, 1:6])
colnames(counts_GV) <- sampleTable$ID[1:6]
counts_GV <- counts_GV[res_MII$sig == "RD" | res_MII$sig == "BLU", ]
chart.Correlation(log(counts_GV + 1))

dds <- DESeqDataSet(se[, 13:18], design = ~Treatment.Control)
colnames(dds) <- sampleTable$ID[13:18]
dds_DESeq <- DESeq(dds)
res_1C <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_1C$sig <- res_MII$sig
counts_1C <- assay(se[, 13:18])
colnames(counts_1C) <- sampleTable$ID[13:18]
counts_1C <- counts_1C[res_MII$sig == "RD" | res_MII$sig == "BLU", ]
chart.Correlation(log(counts_1C + 1))

plot_MA <- function(res, main_lab){
  
  res_plot <- data.frame(mean = res$baseMean, 
                         lfc = res$log2FoldChange, 
                         sig = res$sig)
  
  MA_plot <- ggplot(data = res_plot, aes(x = mean, y = lfc, color = sig, alpha = sig)) + 
    geom_point(size = 0.1) +
    scale_x_log10() +
    scale_y_continuous(limits = c(-4, 4)) + 
    scale_colour_manual(values = c(BLCK = "gray32", RD = "red3", BLU = "blue3")) +
    scale_alpha_manual(values = c(BLCK = 0.3, RD = 1, BLU = 1)) +
    guides(color = FALSE, alpha = FALSE) +
    xlab("mean expression") + 
    ylab("log2FoldChange KO/WT") +
    ggtitle(main_lab)
  return(MA_plot)
}
plot_1C <- plot_MA(res_1C, main_lab = "1C")
plot_MII <- plot_MA(res_MII, main_lab = "MII")
plot_GV <- plot_MA(res_GV, main_lab = "GV")

ggsave("1C_MA_plot.pdf", plot_1C)
ggsave("MII_MA_plot.pdf", plot_MII)
ggsave("GV_MA_plot.pdf", plot_GV)

# fpkm
dds_MII <- DESeqDataSet(se[, 7:12], design = ~1)
fpkm_MII <- fpkm(dds_MII, robust = F)
colnames(fpkm_MII) <- sampleTable$ID[7:12]
write.csv(fpkm_MII, "FPKM/fpkm_MII.csv")

dds_GV <- DESeqDataSet(se[, 1:6], design = ~1)
fpkm_GV <- fpkm(dds_GV, robust = F)
colnames(fpkm_GV) <- sampleTable$ID[1:6]
write.csv(fpkm_GV, "FPKM/fpkm_GV.csv")

dds_1C <- DESeqDataSet(se[, 13:18], design = ~1)
fpkm_1C <- fpkm(dds_1C, robust = F)
colnames(fpkm_1C) <- sampleTable$ID[13:18]
write.csv(fpkm_1C, "FPKM/fpkm_1C.csv")

counts_all <- assay(se) 
colnames(counts_all) <- sampleTable$ID
write.csv(counts_all, "counts_all_samples.csv")

