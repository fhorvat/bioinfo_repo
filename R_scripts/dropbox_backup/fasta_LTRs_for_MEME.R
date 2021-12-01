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
                             mcols(ranges_MT)$repName, "|", 
                             mcols(ranges_MT)$location, "|", 
                             mcols(ranges_MT)$expression), 
              as.string = TRUE,
              file.out = file_name, 
              open = "w")
}

MT_all <- read.delim("MT_MT2_ORR_MLT_properFullLengthLTR.txt", stringsAsFactors = F)

MT <- MT_all
MT <- MT[, c("seqnames", "start", "end", "strand", "repName.x", "rep.location", "s_GV.WE")]
colnames(MT) <- c("seqnames", "start", "end", "strand", "repName", "location", "expression")

MT <- MT[!(grepl("MT2", MT$repName) | grepl("-int", MT$repName)), ]
MT <- MT[grepl("MT", MT$repName), ]

#expressed in oocyte
MT_expressed <- MT[MT$expression > 0, ]
MT_expressed_range <- makeGRangesFromDataFrame(MT_expressed, keep.extra.columns = TRUE)
#getSeqWriteFASTA(MT_expressed_range, "meme/MT_expressed/MT_expressed_in_GV.fasta")

#expressed in oocyte, not MTA/MTB
MT_expressed_notAB <- MT_expressed[!(grepl("MTA", MT_expressed$repName) | grepl("MTB", MT_expressed$repName)), ]
MT_expressed_notAB_range <- makeGRangesFromDataFrame(MT_expressed_notAB, keep.extra.columns = TRUE)
getSeqWriteFASTA(MT_expressed_notAB_range, "meme/MT_expressed/MT_expressed_in_GV_not_AB.fasta")

#expressed in oocyte + 50 not expressed from each class
MT_not_expressed <- MT[MT$expression == 0, ]
MT_not_expressed_sample <- MT_not_expressed %>% 
  group_by(repName) %>%
  sample_n(size = 50)
MT_expressed_plus_50 <- rbind(MT_expressed, MT_not_expressed_sample)
MT_expressed_plus_50_range <- makeGRangesFromDataFrame(MT_expressed_plus_50, keep.extra.columns = TRUE)
getSeqWriteFASTA(MT_expressed_plus_50_range, "meme/MT_expressed/MT_expressed_in_GV_plus_50_from_each_class.fasta")

#expressed in oocyte, not MTA/MTB
MT_expressed_onlyA <- MT_expressed[grepl("MTA", MT_expressed$repName), ]
MT_expressed_onlyA_range <- makeGRangesFromDataFrame(MT_expressed_onlyA, keep.extra.columns = TRUE)
getSeqWriteFASTA(MT_expressed_onlyA_range, "meme/MT_expressed/MT_expressed_in_GV_only_MTA.fasta")
