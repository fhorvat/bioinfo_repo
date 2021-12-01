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
filenames <- file.path(c("STAR/filtered_fasta/Pepa_filter/Fugaku_blastocyst/s_Blast.WE_Aligned.sortedByCoord.out.bam", 
                         "STAR/filtered_fasta/Pepa_filter/Fugaku_GV/s_GV.WE_Aligned.sortedByCoord.out.bam", 
                         "STAR/filtered_fasta/Pepa_filter/Encode_testis/testis_Aligned.sortedByCoord.out.bam"))

bam_Blast <- readGAlignments(filenames[1])
bam_Blast <- subsetByOverlaps(bam_Blast, features)
bam_Blast <- granges(bam_Blast)
se_Blast <- assay(summarizeOverlaps(features = features, reads = bam_Blast, mode = "Union", singleEnd = FALSE, ignore.strand = FALSE))

bam_GV <- readGAlignments(filenames[2])
bam_GV <- subsetByOverlaps(bam_GV, features)
bam_GV <- granges(bam_GV)
se_GV <- assay(summarizeOverlaps(features = features, reads = bam_GV, mode = "Union", singleEnd = FALSE, ignore.strand = FALSE))

bam_testis <- readGAlignments(filenames[3])
bam_testis <- subsetByOverlaps(bam_testis, features)
bam_testis <- granges(bam_testis)
se_testis <- assay(summarizeOverlaps(features = features, reads = bam_testis, mode = "Union", singleEnd = FALSE, ignore.strand = FALSE))

se_all <- data.frame(cbind(se_Blast, se_GV, se_testis))
colnames(se_all) <- c("Fugaku_blastocyst", "Fugaku_GV", "Encode_testis")
rownames(se_all) <- as.character(seqnames(features))

se_all$average <- rowMeans(se_all)
se_all <- se_all[order(se_all$average, decreasing = T), ] 
write.csv(se_all, "STAR/filtered_fasta/Pepa_filter/full_length_L1_Fugaku_Testis_counts.csv")

bed_ordered <- import.bed(con = "references/filtered_fasta/Pepa_filter/full_length_L1_edit.bed")
bed_ordered <- as.data.frame(bed_ordered, row.names = seqnames(features))
bed_ordered <- bed_ordered[rownames(se_all), ]
bed_ordered <- bed_ordered[, c("seqnames", "start", "end", "name", "score", "strand")]
write.table(bed_ordered, "STAR/filtered_fasta/Pepa_filter/full_length_L1_Fugaku_testis_edit_ordered.bed", quote = F, sep = "\t", row.names = F, col.names = F)  
