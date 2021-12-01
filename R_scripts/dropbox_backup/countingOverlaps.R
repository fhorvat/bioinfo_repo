library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)
library(dplyr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/")
features <- import.bed(con = "references/L1_6000bp/L1_6000bp_12column_from_fasta.bed")

filenames <- file.path(paste0("STAR/L1_over_6000bp/mismatch_0.01/11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))
bamfiles <- BamFileList(filenames, yieldSize = 2000000)

register(MulticoreParam())
se <- summarizeOverlaps(features = features, 
                        reads = bamfiles, 
                        mode = "Union", 
                        singleEnd = FALSE, 
                        ignore.strand = FALSE)
colnames(se) <- gsub("_.*", "", colnames(se))
rownames(se) <- paste0(seqnames(features), ":", mcols(features)$name)

se_counts <- as.data.frame(assay(se))
se_counts$average <- rowMeans(se_counts)
se_counts <- se_counts[order(se_counts$average, decreasing = T),] 
# write.csv(se_counts, "STAR/L1_over_6000bp/mismatch_0.2/counts/L1_6000bp_0.2_counts.csv")

dds <- DESeqDataSet(se, design = ~1)
se_fpkm <- as.data.frame(fpkm(dds))
rownames(se_fpkm) <- rownames(se)
se_fpkm$average <- rowMeans(se_fpkm)
se_fpkm <- se_fpkm[order(se_fpkm$average, decreasing = T),] 
# write.csv(se_fpkm, "STAR/L1_over_6000bp/mismatch_0.2/counts/L1_6000bp_0.2_FPKM.csv")