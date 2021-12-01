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
colnames(results_all) <- c("GV", "MII", "oneC")

results_all_one_NA <- results_all[rowSums(is.na(results_all)) < 2, ]
# results_all_no_NA <- results_all[rowSums(is.na(results_all)) < 1, ]

# pheatmap 
# hmcols <- colorRampPalette(c("red", "white", "blue"))(length(bk) - 1)
# bk <- seq(range(results_all_one_NA)[1], range(results_all_one_NA)[2], length = 50)

bk <- unique(c(seq(-2, -0.5, length = 20), seq(-0.5, 0.5, length = 20), seq(0.5, 2, length = 20)))
hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)
pheatmap(results_all_one_NA, 
         clustering_method = "ward.D",
         col = hmcols, 
         breaks = bk,
         fontsize = 15,
         row_clusters = T,
         treeheight_row = 300, 
         show_rownames = F,
         cluster_cols = F)
savepng("pheatmap_all_genes_one_NA_ward.d", width = 1000, asp = 0.8)

# Dynamic Tree Cut
results_all_output_distance <- as.matrix(results_all_one_NA)
results_all_output_distance <- dist(results_all_output_distance)
hc.rows <- hclust(results_all_output_distance, method = "ward.D")
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
savepng("dtc_all_genes_one_NA_ward.d93", width = 1000, asp = 0.75)

results_clusters <- data.frame(results_all_one_NA, dynamicCut)
write.csv(results_all_one_NA_clusters, "dtc_all_genes_one_NA_ward.d93.csv")