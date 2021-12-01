library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)
library(dplyr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/")

features <- import.bed(con = "references/filtered_fasta/Pepa_filter/full_length_L1.bed")
features_names <- paste0(mcols(features)$name, "|", seqnames(features), ":", start(features), "-", end(features), ":", strand(features))
features_end <- width(features)
start(features) <- 1
end(features) <- features_end
features <- as.data.frame(features)
features$seqnames <- features_names
features <- features[, c("seqnames", "start", "end", "name", "score", "strand")]
# write.table(features, "references/filtered_fasta/Pepa_filter/full_length_L1_edit.bed", quote = F, sep = "\t", row.names = F, col.names = F)
features <- makeGRangesFromDataFrame(features)

# bam_files
filenames <- file.path(paste0("STAR/filtered_fasta/Pepa_filter/CNOT6_WT_GV/11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))

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
write.csv(se_all, "STAR/filtered_fasta/Pepa_filter/CNOT6_WT_GV/coverage/full_length_L1_counts.csv")

bed_ordered <- import.bed(con = "references/filtered_fasta/Pepa_filter/full_length_L1_edit.bed")
bed_ordered <- as.data.frame(bed_ordered, row.names = seqnames(features))
bed_ordered <- bed_ordered[rownames(se_all), ]
bed_ordered <- bed_ordered[, c("seqnames", "start", "end", "name", "score", "strand")]
write.table(bed_ordered, "STAR/filtered_fasta/Pepa_filter/CNOT6_WT_GV/coverage/full_length_L1_edit_ordered.bed", quote = F, sep = "\t", row.names = F, col.names = F)  
