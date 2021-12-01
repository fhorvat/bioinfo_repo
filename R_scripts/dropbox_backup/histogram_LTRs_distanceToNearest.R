library("GenomicRanges")
library("Biostrings")
library("BSgenome.Mmusculus.UCSC.mm10")
library("seqinr")
library("dplyr")
library("ggplot2")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/")

elements_all <- read.delim("MT_MT2_ORR_MLT_allElements.txt", header = T, stringsAsFactors = F)
elements <- elements_all
# elements <- elements[!grepl("-int", elements$repName), ]
# elements <- makeGRangesFromDataFrame(elements, keep.extra.columns = TRUE)
# elements_distance <- as.data.frame(distanceToNearest(elements))
# elements_distance$repName <- mcols(elements)$repName[elements_distance$queryHits]

elements_LTR <- elements[!grepl("int", elements$repName), ]
elements_LTR_ranges <- makeGRangesFromDataFrame(elements_LTR, keep.extra.columns = TRUE)
elements_LTR_distance <- as.data.frame(distanceToNearest(elements_LTR_ranges))
elements_LTR_distance$repName <- mcols(elements_LTR_ranges)$repName[elements_LTR_distance$queryHits]
MT2_LTR_distance <- elements_LTR_distance[grepl("MT2", elements_LTR_distance$repName), ]
MT2_LTR_distance <- MT2_LTR_distance[(MT2_LTR_distance$distance) > 99 & (MT2_LTR_distance$distance < 10001), ]

ggplot() + 
  geom_histogram(data = MT2_LTR_distance, aes(distance, fill = repName), binwidth = 20, color = "black") + 
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(repName ~ .)

MT <- elements[grepl("MT2", elements$repName), ]
MT <- MT[!grepl("MT2", MT$repName), ]
MT <- MT[!grepl("-int", MT$repName), ]
MT <- makeGRangesFromDataFrame(MT, keep.extra.columns = TRUE)
# MT <- MT[width(MT) > 300, ]
MT_distance <- as.data.frame(distanceToNearest(MT))
MT_distance$repName <- mcols(MT)$repName[MT_distance$queryHits]

ggplot() + 
  geom_histogram(data = MT2_distance, aes(distance), binwidth = 50, color = "black") + 
  #geom_freqpoly(data = MT2_distance, aes(MT2_distance$distance), binwidth = 20, size = 1.5) +
  scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(repName ~ .)
