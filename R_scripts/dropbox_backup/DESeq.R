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
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("gage")
library("pathview")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/")

writeResults <- function(x, y){
  res_df$gene_name <- mapIds(org.Mm.eg.db,
                            keys = row.names(res_df),
                            column ="GENENAME",
                            keytype = "ENTREZID",
                            multiVals = "first")
  write.csv(res_df, paste0(x, "_vs_", y, "_", sampleTable_exp$Time.Course[1], "_all.csv")) 
  res_df <- subset(res_df, padj < 0.1)
  write.csv(res_df, paste0(x, "_vs_", y, "_", sampleTable_exp$Time.Course[1], "_significant.csv")) 
}

a <- 1:18
filenames <- file.path("/common/WORK/kristian/Projekti/Petr/Cnot6L/Mapping/bbmap/mm9", paste0("11919X", a, "_sorted.bam"))

sampleTable <- read.csv("CNOT6L_sample_list_11919R_2015-10-29.csv", header = T)
sampleTable_exp <- sampleTable[a, ]
# testPairedEndBam(filenames[1])

# KEGG set for mouse
kg.mouse <- kegg.gsets("mouse")
kegg.gs <- kg.mouse$kg.sets[kg.mouse$sigmet.idx]

# counts
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, by = "gene")
register(MulticoreParam())
se <- summarizeOverlaps(features = ebg, 
                        reads = bamfiles, 
                        mode = "Union", 
                        singleEnd = FALSE, 
                        ignore.strand = TRUE)

colData(se) <- DataFrame(sampleTable_exp)

dds <- DESeqDataSet(se, design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID
dds_DESeq <- DESeq(dds)
res <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))

# Plotting 
rld <- rlog(dds, blind = T)
sampleDists <- dist(t(assay(rld)))

# plot distance matrix (heatmap)
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
savepng(paste0(sampleTable_exp$Time.Course[1], "_dist_heatmap_plot"), width = 1000, asp = 0.75)

# PCA
plotPCA(rld, intgroup = c("Treatment.Control", "Time.Course"))  + geom_point(size = 5)
savepng(paste0(sampleTable_exp$Time.Course[1], "_pca_plot"), width = 1000, asp = 0.75)

data <- plotPCA(rld, intgroup = c( "Treatment.Control", "Time.Course"), returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color = Treatment.Control, shape = Time.Course)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
savepng(paste0(sampleTable_exp$Time.Course[1], "_pca_plot_2"), width = 1000, asp = 0.75)

# MA 
plotMA(res, ylim = c(-2, 2), cex = 0.65)
savepng(paste0(sampleTable_exp$Time.Course[1], "_MA_plot"), width = 1000, asp = 0.75)

# write results, make pathview graphs
writeResults("KO", "WT")
gagePathview("KO", "WT")



