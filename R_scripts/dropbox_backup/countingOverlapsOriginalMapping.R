library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)
library(dplyr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/")
features <- import.bed(con = "references/L1_6000bp.bed")
filenames <- file.path(paste0("STAR/original/11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))

bam_x4 <- readGAlignments(filenames[1])
bam_x4 <- subsetByOverlaps(bam_x4, features, minoverlap = 100)
bam_x4 <- granges(bam_x4)
se_x4 <- assay(summarizeOverlaps(features = features, reads = bam_x4, mode = "Union", singleEnd = FALSE, ignore.strand = FALSE))

bam_x5 <- readGAlignments(filenames[2])
bam_x5 <- subsetByOverlaps(bam_x5, features, minoverlap = 100)
bam_x5 <- granges(bam_x5)
se_x5 <- assay(summarizeOverlaps(features = features, reads = bam_x5, mode = "Union", singleEnd = FALSE, ignore.strand = FALSE))

bam_x6 <- readGAlignments(filenames[3])
bam_x6 <- subsetByOverlaps(bam_x6, features, minoverlap = 100)
bam_x6 <- granges(bam_x6)
se_x6 <- assay(summarizeOverlaps(features = features, reads = bam_x6, mode = "Union", singleEnd = FALSE, ignore.strand = FALSE))

se_all <- data.frame(cbind(se_x4, se_x5, se_x6))
colnames(se_all) <- c("11919X4", "11919X5", "11919X6")
rownames(se_all) <- paste(seqnames(features), start(features), end(features), strand(features), sep = ":")

se_all$average <- rowMeans(se_all)
se_all <- se_all[order(se_all$average, decreasing = T), ] 
write.csv(se_all, "STAR/original/L1_6000bp_original_counts.csv")

bed_top10_filtered <- read.delim("references/L1_6000bp.bed", header = F)
colnames(bed_top10_filtered) <- c("seqnames", "start", "end", "name", "score", "strand")
rownames(bed_top10_filtered) <- paste(seqnames(features), start(features), end(features), strand(features), sep = ":")

bed_top10_filtered <- bed_top10_filtered[rownames(se_all), ]
write.table(bed_top10_filtered, "STAR/original/coverage/L1_6000bp_ordered.bed", quote = F, sep = "\t", row.names = F, col.names = F)  
