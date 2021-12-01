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
                             mcols(ranges_MT)$repName.x, "|", 
                             mcols(ranges_MT)$rep.location, "|", 
                             mcols(ranges_MT)$s_GV.WE), 
              as.string = TRUE,
              file.out = file_name, 
              open = "w")
}

elements_vedran <- read.csv("160426.Repeat_Exon_Splicing.csv", stringsAsFactors = F)
elements_vedran <- elements_vedran[, c("rep.coord",  "rep.strand", "repName", "rep.location", "s_GV.WE")]
elements_vedran <- elements_vedran[elements_vedran$s_GV.WE > 0 & !grepl("-int", elements_vedran$repName), ]
elements_vedran$rep.coord <- paste0(elements_vedran$rep.coord, ":", elements_vedran$rep.strand)
elements_vedran <- elements_vedran[, -2]
elements_all <- read.delim("MT_MT2_ORR_MLT_allElements.txt", stringsAsFactors = F)

elements <- elements_all
elements <- elements[!grepl("-int", elements$repName), ]
elements_ranges <- makeGRangesFromDataFrame(elements, keep.extra.columns = TRUE)

width_by <- by(width(elements_ranges), elements_ranges$repName, summary)
width_by <- do.call(rbind, width_by)
width_by <- width_by[, c(2, 5)]

elements_ranges_list <- lapply(X = 1:nrow(width_by), 
                               FUN = function(X) elements_ranges[mcols(elements_ranges)$repName == rownames(width_by)[X]
                                                                 & width(elements_ranges) > width_by[X, 1] 
                                                                 & width(elements_ranges) < width_by[X, 2], ])
elements_ranges_filtered_by_width <- do.call("c", elements_ranges_list)
elements_ranges_filtered_by_width <- sortSeqlevels(elements_ranges_filtered_by_width)
elements_ranges_filtered_by_width <- sort(elements_ranges_filtered_by_width)

elements_ranges_filtered_by_width_df <- as.data.frame(elements_ranges_filtered_by_width)
elements_ranges_filtered_by_width_df$rep.coord <- paste0(elements_ranges_filtered_by_width_df$seqnames, ":", 
                                                         elements_ranges_filtered_by_width_df$start, "-", 
                                                         elements_ranges_filtered_by_width_df$end, ":",
                                                         elements_ranges_filtered_by_width_df$strand)

elements_ranges_filtered_by_width_df_all <- merge(elements_ranges_filtered_by_width_df, elements_vedran, by = "rep.coord", all.x = T, all.y = F)
elements_ranges_filtered_by_width_df_all <- elements_ranges_filtered_by_width_df_all[, c("seqnames", "start", "end", "strand", 
                                                                                         "repName.x", "repName.y", "rep.location", "s_GV.WE")]
elements_ranges_filtered_by_width_df_all$s_GV.WE[is.na(elements_ranges_filtered_by_width_df_all$s_GV.WE)] <- 0
elements_ranges_filtered_by_width_df_all <- makeGRangesFromDataFrame(elements_ranges_filtered_by_width_df_all, keep.extra.columns = TRUE)
getSeqWriteFASTA(elements_ranges_filtered_by_width_df_all, "MT_MT2_ORR_MLT_properFullLengthLTR.fasta")

elements_ranges_filtered_by_width_df_all <- as.data.frame(elements_ranges_filtered_by_width_df_all)
elements_ranges_filtered_by_width_df_all <- elements_ranges_filtered_by_width_df_all[, c("seqnames", "start", "end", "strand", "repName.x", "rep.location", "s_GV.WE")]
write.table(elements_ranges_filtered_by_width_df_all, "MT_MT2_ORR_MLT_properFullLengthLTR.txt", quote = F, sep = "\t", row.names = F, col.names = T)
