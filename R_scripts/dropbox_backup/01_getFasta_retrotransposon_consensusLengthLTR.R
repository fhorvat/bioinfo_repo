library("GenomicRanges")
library("Biostrings")
library("BSgenome.Mmusculus.UCSC.mm10")
library("seqinr")
library("dplyr")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/")

getSeqWriteFASTA <- function(ranges_MT, file_name){
  ranges_MT_seq <- getSeq(Mmusculus,
                          seqnames(ranges_MT), 
                          start(ranges_MT), 
                          end(ranges_MT))
  ranges_MT_seq <- as.character(ranges_MT_seq)
  write.fasta(as.list(ranges_MT_seq),
              nbchar = 80, 
              names = paste0(seqnames(ranges_MT), ":", 
                             start(ranges_MT), "-", 
                             end(ranges_MT), ":",
                             strand(ranges_MT), "|", 
                             mcols(ranges_MT)$repName), 
              as.string = TRUE,
              file.out = file_name, 
              open = "w")
}

# getting retrotransposon classes and consensus lengths
elements_length <- read.csv("../retrotransposon_conservation/retrotransposon_families_characteristics.csv", stringsAsFactors = F)
elements_length <- elements_length[elements_length$Type != 'int', ]
elements_length <- elements_length[, c(1, 2, 8)]
colnames(elements_length) <- c("retrotransposonID", "LTRClass", "consensusLength")

# getting all retrotransposon elements coordinates, filtering -int sequences
elements_all <- read.delim("MT_MT2_ORR_MLT_allElements.txt", stringsAsFactors = F)
elements_all <- elements_all[!grepl("-int", elements_all$repName), ]

# merging all LTRs with classes/consensus lengths, adding coordinates as name for merging with Vedran's table
elements_all_merged <- merge(x = elements_all, y = elements_length, by.x = "repName", by.y = "retrotransposonID")
elements_all_merged <- elements_all_merged[, c(2:5, 1, 6, 7)]
elements_all_merged$repCoord <- paste0(elements_all_merged$seqnames, ":", 
                                       elements_all_merged$start, "-", 
                                       elements_all_merged$end, ":",
                                       elements_all_merged$strand)

# filtering by consesus length
elements_all_merged_filtered <- elements_all_merged
elements_all_merged_filtered$width <- elements_all_merged$end - elements_all_merged$start + 1
elements_all_merged_filtered <- elements_all_merged_filtered[, c(1:6, 8, 7, 9)]

# # +-50 nt 
# elements_all_merged_GV_filtered$in_length <- ifelse((elements_all_merged_filtered$width >= (elements_all_merged_filtered$consensusLength - 50)) &
#                                                       (elements_all_merged_filtered$width <= (elements_all_merged_filtered$consensusLength + 50)), 
#                                                     T, F)

# +-5% length 
elements_all_merged_filtered$length_5_percent <- elements_all_merged_filtered$consensusLength * 0.05
elements_all_merged_filtered$in_length <- ifelse((elements_all_merged_filtered$width >= 
                                                       (elements_all_merged_filtered$consensusLength - elements_all_merged_filtered$length_5_percent)) &
                                                      (elements_all_merged_filtered$width <= 
                                                         (elements_all_merged_filtered$consensusLength + elements_all_merged_filtered$length_5_percent)), 
                                                    T, F)

elements_all_merged_filtered <- elements_all_merged_filtered[elements_all_merged_filtered$in_length == T, ]
elements_all_merged_filtered_GR <- makeGRangesFromDataFrame(elements_all_merged_filtered, keep.extra.columns = TRUE)
getSeqWriteFASTA(elements_all_merged_filtered_GR, "../retrotransposon_conservation/retrotransposon_consensusLength_5percent_LTR.fasta")


