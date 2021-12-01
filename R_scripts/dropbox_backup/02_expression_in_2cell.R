library("GenomicRanges")
library("Biostrings")
library("seqinr")
library("dplyr")
library("GenomicAlignments")
library("dplyr")
library("DataCombine")
library("DESeq2")
library("tidyr")
library("rtracklayer")
library("CoverageView")
library("IRanges")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/coverage")

# reading 
MT2_solo <- read.delim("MT2_solo_LTRs_Maja_20160803.txt", stringsAsFactors = F)
ORR1A0_solo <- read.delim("ORR1A0_solo_LTRs_Maja_20160803.txt", stringsAsFactors = F)

# combining both repeats, making GRanges
MT2_ORR1A0_solo <- rbind(MT2_solo, ORR1A0_solo)
MT2_ORR1A0_solo_gr <- makeGRangesFromDataFrame(MT2_ORR1A0_solo, keep.extra.columns = T)

# making TxDb object from knownGene gtf from UCSC
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mm10/UCSC_knownGenes_mm10_20160513.gtf.gz")
knownGenes_gtf_gr <- genes(knownGenes_gtf)

# finding overlaps between MT2/ORR1A0 solo LTRs and knownGene table, filtering those which overlap
MT2_ORR1A0_solo_knownGenes_overlaps <- findOverlaps(MT2_ORR1A0_solo_gr, knownGenes_gtf_gr)
MT2_ORR1A0_solo_filtered <- MT2_ORR1A0_solo_gr[-queryHits(MT2_ORR1A0_solo_knownGenes_overlaps), ]

#### calculating FPKMs of MT2/ORR1A0 solo LTRs
# bam filepaths
bam_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))

# library size
logs_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out", 
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAmLog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WELog.final.out"))

# calculating library size in millions of reads
number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- c("GV", "1C", "2C", "2C_aphi", "4C")
number_of_reads <- number_of_reads / 10^6

# counts for 2 cell stage
bamfiles <- BamFileList(bam_filenames, yieldSize = 2000000)
se <- summarizeOverlaps(features = MT2_ORR1A0_solo_filtered, 
                        reads = bamfiles["s_2cell.WE.bam"], 
                        mode = "Union", 
                        singleEnd = FALSE, 
                        ignore.strand = TRUE)

# data.frame with counts 
count_df <- as.data.frame(assay(se))
colnames(count_df) <- "bam_2C"

# calculating fpkm from counts
fpkm_df <- count_df
fpkm_df$width <- width(MT2_ORR1A0_solo_filtered)
MT2_ORR1A0_solo_filtered$FPKM <- fpkm_df$bam_2C / (number_of_reads["2C"] * (fpkm_df$width / 1000))

# order by FPKM in 2cell
MT2_ORR1A0_solo_filtered_df <- as.data.frame(MT2_ORR1A0_solo_filtered)
MT2_ORR1A0_solo_filtered_df_ordered <- 
  MT2_ORR1A0_solo_filtered_df %>%
  group_by(repName) %>%
  arrange(desc(FPKM))

# writting ordered table
write.table(MT2_ORR1A0_solo_filtered_df_ordered, "MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", quote = F, sep = "\t", row.names = F, col.names = T)


