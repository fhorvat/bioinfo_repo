# ###############################
# ### 4 BEFORE WE START
# ###############################
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2", "genefilter")
# install.packages("pheatmap")

library("airway")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")
library("AnnotationDbi")
library("org.Hs.eg.db")

##############################
### 6 LOCATING ALIGMENT FILES
###############################
dir <- system.file("extdata", package = "airway", mustWork = TRUE)
# dir <- "/home/manager/R/x86_64-pc-linux-gnu-library/3.2/airway/extdata"
list.files(dir)
bam_filenames <- list.files(dir, pattern = '.bam', full.names = T)
sample_table <- read.csv(paste0(dir, '/sample_table.csv'), row.names = 1)


###############################
### 7 READING .BAM FILES
###############################
bamfiles <- BamFileList(bam_filenames, yieldSize = 2000000)
seqinfo(bamfiles[1]) 


###############################
### 8 READING .BAM FILES
###############################
gtffile <- paste0(dir, "/Homo_sapiens.GRCh37.75_subset.gtf")

# making TxDb object from GTF
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
asGFF(txdb)


###############################
### 9 GROUPING EXONS BY GENES
###############################
ebg <- exonsBy(txdb, by = "gene")
head(as.data.frame(ebg))


###############################
### 10 COUNTING READS
###############################
se <- summarizeOverlaps(features = ebg, 
                        reads = bamfiles,
                        mode = "Union",
                        singleEnd = FALSE,
                        ignore.strand = TRUE,
                        fragments = TRUE) # fragments = T counts unpaired reads in paired end experiment
se
dim(se)
class(se)


###############################
### 12 COUNTING READS
###############################
# assay (pink block) - matrix of counts 
head(assay(se))

# rowRanges (blue block) - information about the genomic ranges 
rowRanges(se)
str(metadata(rowRanges(se)))

# colData (green block) - information about the samples
colData(se)
(colData(se) <- DataFrame(sample_table))

se$cell
se$dex 


###############################
### 13 EXPLORING COUNTS
###############################
# reading in sample SummarisedExperiment object 
data("airway")
se <- airway
head(assay(se))

# quick check the millions of fragments that uniquely aligned to the genes
round(colSums(assay(se)) / 1e6, 1)
colData(se)

# quick info
str(metadata(rowRanges(se)))
summary(assay(se))


###############################
### 14 STARTING DESeq2
###############################
## DESeqDataSet
dds <- DESeqDataSet(se, design = ~ cell + dex)
head(counts(dds))


###############################
### 15 DESeq2
###############################
# pre-filtering the dataset
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)

# getting fpkm 
dds_fpkm <- fpkm(dds)
head(dds_fpkm)


###############################
### 17 EXPLORATORY ANALYSIS AND VISUALIZATION
###############################
#  rlog transformation
rld <- rlog(dds, blind = F)
head(assay(rld))

# plotting difference between log2 + 1 and rlog
par(mfrow = c(1, 2))
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized = TRUE)[,1:2] + 1), pch = 16, cex = 0.3)
plot(assay(rld)[,1:2], pch = 16, cex = 0.3)


###############################
### 18 EXPLORATORY ANALYSIS AND VISUALIZATION
###############################
# sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists


###############################
### 19 EXPLORATORY ANALYSIS AND VISUALIZATION
###############################
# sample distances visualization
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$dex, rld$cell, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


###############################
### 20 EXPLORATORY ANALYSIS AND VISUALIZATION
###############################
# PCA plot
plotPCA(rld, intgroup = "dex")
plotPCA(rld, intgroup = "cell")


###############################
### 21 EXPLORATORY ANALYSIS AND VISUALIZATION
###############################
# PCA plot ggplot2
data <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


###############################
### 22 DIFFERENTIAL EXPRESSION ANALYSIS
###############################
# Running the differential expression pipeline
dds <- DESeq(dds)
# building the results table
(res <- results(dds))


###############################
### 23 DIFFERENTIAL EXPRESSION ANALYSIS
###############################
# metadata of results table (or what each column represents?)
mcols(res, use.names = TRUE)


###############################
### 25 EXPLORING RESULTS
###############################
# summary
summary(res)

# p-adjusted = 0.05
res.05 <- results(dds, alpha = .05)
table(res.05$padj < 0.05)

# lfcThreshold = 1
resLFC1 <- results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.1)


###############################
### 26 EXPLORING RESULTS
###############################
# getting different relations in results
results(dds, contrast = c("cell", "N061011", "N61311"))


###############################
### 27 EXPLORING RESULTS
###############################
# strongest downregulation
resSig <- subset(res, padj < 0.1)
head(resSig[order(resSig$log2FoldChange), ])


###############################
### 28 EXPLORING RESULTS
###############################
# strongest upregulation
head(resSig[order(resSig$log2FoldChange, decreasing = T), ])


###############################
### 29 PLOTTING RESULTS
###############################
# visualize the counts for a particular gene
topGene <- rownames(res)[which.min(res$padj)]
dev.off()
plotCounts(dds, gene = topGene, intgroup = c("dex"), pch = 20)


###############################
### 30 PLOTTING RESULTS
###############################
# visualize the counts for a particular gene - ggplot2
data <- plotCounts(dds, gene = topGene, intgroup = c("dex", "cell"), returnData = TRUE)
ggplot(data, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() + 
  geom_point(position = position_jitter(width = 0.1, height = 0), size = 3)


###############################
### 31 PLOTTING RESULTS
###############################
# maplot
plotMA(res, ylim = c(-5,5))


###############################
### 32 PLOTTING RESULTS
###############################
# gene clustering
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = T), 20)
mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[, c("cell", "dex")])
pheatmap(mat, annotation_col = df)


###############################
### 33 ANNOTATING AND EXPORTING RESULTS
###############################
# available key columns
columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$entrez <- mapIds(org.Hs.eg.db, keys = row.names(res), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
res$gene_name <- mapIds(org.Hs.eg.db, keys = row.names(res), column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
head(res)


###############################
### 34 ANNOTATING AND EXPORTING RESULTS
###############################
#exporting results
res_ordered <- res[order(res$log2FoldChange), ]
res_ordered <- as.data.frame(res_ordered)[1:100,]
write.csv(res_ordered, file = "results.csv")
getwd()
