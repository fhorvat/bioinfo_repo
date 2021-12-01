library("GenomicAlignments")
library("Rsamtools")
library("BiocParallel")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/Fugaku")

data_path <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2"
track_filenames <- list.files(path = data_path, 
                              pattern = "*.bam$", 
                              recursive = T, 
                              full.names = T)

logs_filenames <- list.files(path = data_path, 
                             pattern = "*Log.final.out", 
                             recursive = T, 
                             full.names = T)

sample_names <- gsub("^/.*/|\\.bam", "", track_filenames)

# Uniquely mapped reads number in millions
number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- sample_names
number_of_reads <- number_of_reads / 10^6

# counts over exons of all genes from knownGene table from UCSC
bamfiles <- BamFileList(track_filenames, yieldSize = 2000000)
ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
register(MulticoreParam())
se_Fugaku <- summarizeOverlaps(features = ebg, 
                               reads = bamfiles, 
                               mode = "Union", 
                               singleEnd = FALSE, 
                               ignore.strand = TRUE)

count_df <- as.data.frame(assay(se_Fugaku))
colnames(count_df) <- sample_names

# calculating fpkm from counts  
fpkm_df <- count_df
fpkm_df$width <- sapply(width(ebg), sum)

invisible(lapply(X = names(number_of_reads),
                 FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (number_of_reads[X] * (fpkm_df$width / 1000))))
write.csv(fpkm_df, "fpkm_fugaku_all_samples.csv", quote = F, row.names = T)
write.csv(count_df, "counts_fugaku_all_samples.csv", quote = F, row.names = T)

