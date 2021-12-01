# GO entrichment on all genes in all 3 timepoints (log2FC values KO vs. WT after DESeq2 analysis)
# clustered with DynamicTreeCuter by ward.D method in 15 clusters and cut 

# library("DESeq2")
# library("Rsamtools")
# library("GenomicFeatures")
# library("GenomicAlignments")
# library("BiocParallel")
# library("geneplotter")
# library("pheatmap")
# library("RColorBrewer")
# library("ggplot2")
# library("genefilter")
# library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("dynamicTreeCut")
library("WGCNA")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("topGO")
library("plyr")
library("GOstats")
library("goseq")

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

# dds <- DESeqDataSet(se[, 1:6], design = ~Treatment.Control)
# dds <- dds[rowSums(counts(dds)) > 1, ]
# colnames(dds) <- sampleTable_exp$ID[1:6]
# dds_DESeq <- DESeq(dds)
# res_GV <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
# res_GV <- data.frame(res_GV)
# 
# dds <- DESeqDataSet(se[, 7:12], design = ~Treatment.Control)
# dds <- dds[rowSums(counts(dds)) > 1, ]
# colnames(dds) <- sampleTable_exp$ID[7:12]
# dds_DESeq <- DESeq(dds)
# res_MII <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
# res_MII <- data.frame(res_MII)

# dds <- DESeqDataSet(se[, 13:18], design = ~Treatment.Control)
# dds <- dds[rowSums(counts(dds)) > 1, ]
# colnames(dds) <- sampleTable_exp$ID[13:18]
# dds_DESeq <- DESeq(dds)
# res_1C <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
# res_1C <- data.frame(res_1C)

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

results_all_output_distance <- as.matrix(res_all_significant)
results_all_output_distance <- dist(results_all_output_distance)
hc.rows <- hclust(results_all_output_distance, method = "ward.D")

# Dynamic Tree Cut
maxCoreScatter <- 0.95
minGap <- (1 - maxCoreScatter) * 3/4
dynamicCut <- cutreeDynamic(hc.rows,
                            method = "hybrid", 
                            distM = as.matrix(results_all_output_distance), 
                            deepSplit = 2, 
                            maxCoreScatter = maxCoreScatter, 
                            minGap = minGap, 
                            maxAbsCoreScatter = NULL, 
                            minAbsGap = NULL)
# cut2colour <- labels2colors(dynamicCut) 
# plotDendroAndColors(hc.rows, 
#                     cut2colour, 
#                     "Dynamic Tree Cut", 
#                     dendroLabels = FALSE,
#                     hang = 0.03,
#                     addGuide = TRUE,
#                     guideHang = 0.05) 
# savepng("dtc_clusters_significant_genes_log2FC_wardD", width = 1000)

# significant_genes_log2FC_clusters <- data.frame(res_all_significant, cluster = dynamicCut)
# write.csv(significant_genes_log2FC_clusters, "significant_genes_log2FC_clusters_wardD.csv")

## pheatmap with cluster annotation
# ann_row <- data.frame(cluster = factor(dynamicCut), row.names = rownames(res_all_significant))
# cluster_colors <- labels2colors(1:40) 
# names(cluster_colors) <- 1:40 
# cluster_colors <- cluster_colors[1:length(unique(cut2colour))]
# ann_colors <- list(cluster = cluster_colors)

# bk <- unique(c(seq(-2, -0.2, length = 10), seq(-0.2, 0.2, length = 20), seq(0.2, 2, length = 10)))
# hmcols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(bk) - 1)

# pheatmap_output <- pheatmap(res_all_significant, 
#                             clustering_method = "ward.D",
#                             col = hmcols, 
#                             breaks = bk,
#                             fontsize = 15,
#                             annotation_row = ann_row, 
#                             annotation_colors = ann_colors,
#                             row_clusters = T,
#                             treeheight_row = 300, 
#                             show_rownames = F,
#                             cluster_cols = F)
# savepng("heatmap_clusters_significant_genes_log2FC_wardD", width = 1000)


# GOstats - ENTREZ ID as rownames
GOstats_function <- function(results_clusters, x){
  selectedIDs <- rownames(subset(results_clusters, cluster == x))
  
  goResults <- function(ontology){
    goParams <- new("GOHyperGParams",
                    geneIds = selectedIDs,
                    universeGeneIds = universeIDs,
                    annotation = "org.Mm.eg",
                    ontology = ontology,
                    pvalueCutoff = 0.01,
                    conditional = TRUE,
                    testDirection = "over")
    goResults <- hyperGTest(goParams)
    GOstats_results <- summary(goResults)
    colnames(GOstats_results) <- c("GO.ID", colnames(GOstats_results)[2:7])
    return(GOstats_results)
  }
  
  GOstats_results <- lapply(c("MF", "BP", "CC"), function(y) goResults(y))
  GOstats_results <- rbind.fill(GOstats_results)
  GOstats_results <- GOstats_results[, c(1, 7, 2:6)]
  GOstats_results <- GOstats_results[order(GOstats_results$Pvalue), ]
  return(GOstats_results)
}

results_clusters <- data.frame(res_all_significant, dynamicCut)
results_clusters <- data.frame(results_clusters$dynamicCut, row.names = rownames(results_clusters))
colnames(results_clusters) <- "cluster"
universeIDs <- rownames(results_all_no_NA)

GOstats_results <- lapply(1:max(results_clusters$cluster), function(x) GOstats_function(results_clusters, x))
lapply(1:length(GOstats_results), function(i) write.csv(GOstats_results[[i]], 
                                                        file = paste0("GOstatsResults_cluster", i, "_significant_genes.csv"),
                                                        row.names = FALSE))

