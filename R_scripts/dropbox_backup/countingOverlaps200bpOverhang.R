library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)
library(dplyr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/")
features <- import.bed(con = "references/filtered_fasta/200bp_overhang/L1_6000bp_12column_from_fasta.bed")

# removing 200bp flanking regions from bed features
start(features) <- start(features) + 200
end(features) <- end(features) - 200

# filtering bam files based on overlaping with features from bed
# min. overlap between each read and feature must be at least 100 bp
filenames <- file.path(paste0("STAR/filtered_fasta/200bp_overhang/mismatch_0.01/11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))

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
rownames(se_all) <- as.character(seqnames(features))

se_all$average <- rowMeans(se_all)
se_all <- se_all[order(se_all$average, decreasing = T), ] 
write.csv(se_all, "STAR/L1_over_6000bp/mismatch_0.01/coverage/filtered/L1_6000bp_0.01_counts.csv")

bed_top10_filtered <- read.delim("references/L1_6000bp/L1_6000bp_12column_from_fasta.bed", header = F)
bed_top10_filtered <- bed_top10_filtered[, 1:6]
colnames(bed_top10_filtered) <- c("seqnames", "start", "end", "name", "score", "strand")
rownames(bed_top10_filtered) <- bed_top10_filtered$seqnames
bed_top10_filtered <- bed_top10_filtered[rownames(se_all), ]
write.table(bed_top10_filtered, "STAR/L1_over_6000bp/mismatch_0.01/coverage/filtered/L1_6000bp_ordered.bed", quote = F, sep = "\t", row.names = F, col.names = F)  

bed_top10_filtered <- bed_top10_filtered[rownames(se_all)[1:10], ]
write.table(bed_top10_filtered, "STAR/L1_over_6000bp/mismatch_0.01/coverage/filtered/L1_6000bp_top10.bed", quote = F, sep = "\t", row.names = F, col.names = F)  
