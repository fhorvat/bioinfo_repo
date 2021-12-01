library("dynamicTreeCut")
library("WGCNA")
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
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

# a <- 1:18
# filenames <- file.path("/common/WORK/kristian/Projekti/Petr/Cnot6L/Mapping", paste0("11919X", a, "_sorted.bam"))
# sampleTable <- read.csv("/common/WORK/fhorvat/RNA_Seq_CNOT6L/CNOT6L_sample list_11919R_2015-10-29.csv", header = T)
# sampleTable_exp <- sampleTable[a, ]
# bamfiles <- BamFileList(filenames, yieldSize = 2000000)
# ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
# register(MulticoreParam())
# se <- summarizeOverlaps(features = ebg, 
#                         reads = bamfiles, 
#                         mode = "Union", 
#                         singleEnd = TRUE, 
#                         ignore.strand = TRUE)
# colData(se) <- DataFrame(sampleTable_exp)

dds <- DESeqDataSet(se[, 1:6], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[1:6]
dds_DESeq <- DESeq(dds)
res_GV <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_GV <- data.frame(res_GV)

dds <- DESeqDataSet(se[, 7:12], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[7:12]
dds_DESeq <- DESeq(dds)
res_MII <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_MII <- data.frame(res_MII)

dds <- DESeqDataSet(se[, 13:18], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[13:18]
dds_DESeq <- DESeq(dds)
res_1C <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_1C <- data.frame(res_1C)
                     
results_all <- merge(data.frame(res_GV$log2FoldChange, row.names = rownames(res_GV)),
                     data.frame(res_MII$log2FoldChange, row.names = rownames(res_MII)),
                     by = 0, all = T)
rownames(results_all) <- results_all$Row.names
results_all <- results_all[, -1]
colnames(results_all) <- c("GV", "MII")

results_all <- merge(results_all, 
                     data.frame(res_1C$log2FoldChange, row.names = rownames(res_1C)),
                     by = 0, all = T)
rownames(results_all) <- results_all$Row.names
results_all <- results_all[, -1]
colnames(results_all) <- c("GV", "MII", "1C")

results_all_no_NA <- results_all[rowSums(is.na(results_all)) < 1, ]

# pheatmap 
bk <- unique(c(seq(-2, -0.2, length = 10), seq(-0.2, 0.2, length = 30), seq(0.2, 2, length = 10)))
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)

pheatmap_method <- function(results_all_no_NA, method){
  pheatmap(results_all_no_NA,
           clustering_method = method,
           main = paste("log2FC DESeq2 all genes", method, "method"),
           col = hmcols, 
           breaks = bk,
           cluster_row = T, 
           cluster_cols = F,
           treeheight_row = 100, 
           treeheight_col	= 100,
           show_rownames = F)
  savepng(paste0("pheatmap_log2FC_all_genes_", method), width = 1000)  
}
lapply(c("ward.D", "ward.D2", "complete", "average"), pheatmap_method, results_all_no_NA = results_all_no_NA)

