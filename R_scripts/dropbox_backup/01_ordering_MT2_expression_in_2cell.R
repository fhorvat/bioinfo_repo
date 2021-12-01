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

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/")

load("gr_mis.Robj")
maja_table <- as.data.frame(gr_mis)[, c("seqnames", "start", "end", "strand", "repName", "id")]
maja_table <- maja_table[grepl("MT2", maja_table$repName), ]
maja_table <- maja_table[maja_table$id != 0, ]
maja_ranges_by <- maja_table
maja_ranges_by <- by(maja_ranges_by[, c("start", "end")], maja_ranges_by$id, range)
maja_width <- unlist(lapply(maja_ranges_by, function(x) x[2] - x[1]))
maja_width <- maja_width[maja_width > 5900 & maja_width < 6600]
maja_table_filtered <- maja_table[maja_table$id %in% names(maja_width), ]
maja_gr <- makeGRangesFromDataFrame(maja_table_filtered, keep.extra.columns = TRUE)
MT2_ranges_filtered <- sortSeqlevels(maja_gr)
MT2_ranges_filtered <- sort(MT2_ranges_filtered)
MT2_df_filtered <- as.data.frame(MT2_ranges_filtered)
MT2_df_filtered <- slide(MT2_df_filtered, Var = "end", slideBy = 1)
MT2_df_filtered <- slide(MT2_df_filtered, Var = "id", slideBy = 1)
MT2_df_filtered$distance <- MT2_df_filtered$end1 - MT2_df_filtered$start
MT2_df_filtered <- MT2_df_filtered[MT2_df_filtered$distance < 7000 & MT2_df_filtered$distance > 0, ]
MT2_df_filtered <- head(MT2_df_filtered, -1)
MT2_df_filtered <- MT2_df_filtered[MT2_df_filtered$id == MT2_df_filtered$id1, ]
MT2_df_filtered <- MT2_df_filtered[, c("seqnames", "start", "end1", "strand", "repName", "id")]
colnames(MT2_df_filtered) <- c("seqnames", "start", "end", "strand", "repName", "id")

MT2_ranges_full <- makeGRangesFromDataFrame(MT2_df_filtered, keep.extra.columns = TRUE)
MT2_ranges_full <- MT2_ranges_full[width(MT2_ranges_full) > 6000, ]
# write.table(as.data.frame(MT2_ranges_full), "MT2_fullLengthElements.txt", quote = F, sep = "\t", row.names = F, col.names = F)

MT2_ranges_LTR <- makeGRangesFromDataFrame(maja_table, keep.extra.columns = TRUE)
MT2_ranges_LTR <- MT2_ranges_LTR[mcols(MT2_ranges_LTR)$id %in% mcols(MT2_ranges_full)$id, ]
# write.table(as.data.frame(MT2_ranges_LTR), "MT2_LTRFromFullLengthElements.txt", quote = F, sep = "\t", row.names = F, col.names = F)

bam_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                             "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))

### calculating FPKMs of LTRs
# library size
logs_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out", 
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAmLog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WELog.final.out"))

number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- c("GV", "1C", "2C", "2C_aphi", "4C")
number_of_reads <- number_of_reads / 10^6

# counts 
bamfiles <- BamFileList(bam_filenames, yieldSize = 2000000)
features_filtered <- MT2_ranges_LTR
se <- summarizeOverlaps(features = features_filtered, 
                        reads = bamfiles["s_2cell.WE.bam"], 
                        mode = "Union", 
                        singleEnd = FALSE, 
                        ignore.strand = TRUE)

count_df <- as.data.frame(assay(se))
colnames(count_df) <- "bam_2C"

# calculating fpkm from counts
MT2_LTR_fpkm <- count_df
MT2_LTR_fpkm$width <- width(MT2_ranges_LTR)
MT2_LTR_fpkm[, "bam_2C"] <- MT2_LTR_fpkm$bam_2C / (number_of_reads["2C"] * (MT2_LTR_fpkm$width / 1000))

# ordering by fpkm in 2cell, getting unique LTR IDs
MT2_LTR_fpkm$id <- mcols(features_filtered)$id
MT2_LTR_fpkm_ordered <- MT2_LTR_fpkm[order(MT2_LTR_fpkm[, "bam_2C"], decreasing = T), ]
MT2_LTR_fpkm_ordered_id <- unique(MT2_LTR_fpkm_ordered$id)

# getting ranges of full MT2 ordered by expression of LTRs
MT2_ranges_full_ordered <- MT2_ranges_full[mcols(MT2_ranges_full)$id %in% MT2_LTR_fpkm_ordered_id, ]
MT2_ranges_full_ordered <- MT2_ranges_full_ordered[match(MT2_LTR_fpkm_ordered_id, mcols(MT2_ranges_full_ordered)$id), ] 

write.table(as.data.frame(MT2_ranges_full_ordered), "MT2_fullLengthOrderedByExpressionIn2cell_new.txt", quote = F, sep = "\t", row.names = F, col.names = F)


