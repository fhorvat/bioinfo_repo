library("dplyr")
library("tidyr")
library("GenomicRanges")
library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("BSgenome.Mmusculus.UCSC.mm10")
library("seqinr")
library("DECIPHER")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/motif_discovery/MT/")

# get sequence from ranges and write .fasta
function(ranges_MT, file_name){
  
  ranges_MT_seq <- getSeq(x = Mmusculus,
                          names = ranges_MT$chr,
                          start = as.integer(ranges_MT$start),
                          end = as.integer(ranges_MT$end),
                          strand = ranges_MT$rep.strand)
  
  ranges_MT_seq_oriented <- OrientNucleotides(myXStringSet = ranges_MT_seq,
                                              reference = which.max(width(ranges_MT_seq)),
                                              type = "sequences",
                                              orientation = "all",
                                              threshold = 0.05,
                                              verbose = TRUE,
                                              processors = 4)
  
  ranges_MT_seq <- as.character(ranges_MT_seq)
  write.fasta(as.list(ranges_MT_seq),
              nbchar = 80,
              names = paste0(ranges_MT$chr, ":",
                             ranges_MT$start, "-",
                             ranges_MT$end, ":",
                             ranges_MT$rep.strand, "|",
                             ranges_MT$repName),
              as.string = TRUE,
              file.out = file_name,
              open = "w")
  
  return(file_name)
}


# get sequence from ranges, divide by expression to top/mid/bottom and write .fasta
getSeqWriteFASTAquantiles <- function(ranges_MT, file_name){
  
#   ranges_MT <- selected_MT_fpkm[selected_MT_fpkm$repName == "MTA", ]
 
   ranges_MT <- ranges_MT %>%
    mutate(quantile = ntile(s_GV.WE, 3)) %>%
#     mutate(quantile = replace(quantile, quantile != 1, 0)) %>%
    mutate(quantile = factor(quantile, labels = c("bottom_33", "mid_33", "top_33")))

    invisible(sapply(X = unique(ranges_MT$quantile), 
                   FUN = function(X){
                     
                     # filter by expression quantile
                     ranges_MT <- ranges_MT %>%
                       filter(quantile == X)
                     
                     # get sequence
                     ranges_MT_seq <- getSeq(x = Mmusculus,
                                             names = ranges_MT$chr, 
                                             start = as.integer(ranges_MT$start),  
                                             end = as.integer(ranges_MT$end), 
                                             strand = ranges_MT$rep.strand) 
                     
                     # write.fasta
                     ranges_MT_seq <- as.character(ranges_MT_seq)
                     write.fasta(as.list(ranges_MT_seq),
                                 nbchar = 80, 
                                 names = paste0(ranges_MT$chr, ":", 
                                                ranges_MT$start, "-", 
                                                ranges_MT$end, ":",
                                                ranges_MT$rep.strand, "|", 
                                                ranges_MT$repName, "| FPKM = ", 
                                                ranges_MT$s_GV.WE),  
                                 as.string = TRUE,
                                 file.out = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/motif_discovery/MT/sequences/", 
                                                   file_name, "/", file_name, "_", X, "_fpkm_GV.fasta"),
                                 open = "w")
           }))        
  
  return(file_name)
}

# full contribution to a 5’ exon in protein-coding genes and lncRNAs
repeats_list <- readRDS("/common/WORK/vfranke/Projects/PSvoboda_MT/Results/MT_FindRepeats/mm/160907.mm.Results.rds", refhook = NULL)
selected_MT <- repeats_list$RepeatsSelected_MALR
selected_MT <- 
  selected_MT %>% 
  filter(grepl("MT", repName) & !grepl("MT2|-int", repName)) %>%
  filter(grepl("5' exon", ex.category)) %>%
  filter(grepl("Complete", ex.category.complete))

# select column, separate rep.coord to new columns
selected_MT_df <-  selected_MT %>%
  dplyr:::select(repName, exon_id, rep.coord, rep.strand) %>%
  separate(rep.coord, c("chr", "start", "end"), ":|-", remove = F)

# make GRanges for count features
selected_MT_gr <- makeGRangesFromDataFrame(selected_MT_df, 
                                           start.field = "start",
                                           end.field = "end", 
                                           seqnames.field = "chr",
                                           strand.field = "rep.strand",
                                           keep.extra.columns = T)

# .bam files paths
data_path <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/"
track_filenames <- file.path(paste0(data_path, c("s_1cell.PA/s_1cell.PA.bam",
                                                 "s_MII.PA/s_MII.PA.bam", 
                                                 "s_1cell.WE/s_1cell.WE.bam", 
                                                 "s_GV.WE/s_GV.WE.bam", 
                                                 "s_MII.WE/s_MII.WE.bam")))

# counting Fugaku's data over MT as features 
bamfiles <- BamFileList(track_filenames, yieldSize = 2000000)
register(MulticoreParam())
se_Fugaku <- summarizeOverlaps(features = selected_MT_gr, 
                               reads = bamfiles, 
                               mode = "Union", 
                               singleEnd = FALSE, 
                               ignore.strand = TRUE)

# count data.frame
count_df <- as.data.frame(assay(se_Fugaku))
colnames(count_df) <- c("s_1cell.PA", "s_MII.PA", "s_1cell.WE", "s_GV.WE", "s_MII.WE")
count_df$rep.coord <- selected_MT_gr$rep.coord

### calculating fpkm from counts  
# library size
logs_filenames <- file.path(paste0(data_path, c("s_1cell.PA/s_1cell.PALog.final.out",
                                                "s_MII.PA/s_MII.PALog.final.out",
                                                "s_1cell.WE/s_1cell.WELog.final.out", 
                                                "s_GV.WE/s_GV.WELog.final.out",
                                                "s_MII.WE/s_MII.WELog.final.out")))

# Uniquely mapped reads number in millions
number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- c("s_1cell.PA", "s_MII.PA", "s_1cell.WE", "s_GV.WE", "s_MII.WE")
number_of_reads <- number_of_reads / 10^6

# FPKM data frame
fpkm_df <- count_df
fpkm_df$width <- width(selected_MT_gr)
invisible(lapply(X = names(number_of_reads),
                 FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (number_of_reads[X] * (fpkm_df$width / 1000))))

# join FPKM with coordinate data
selected_MT_fpkm <- left_join(selected_MT_df, fpkm_df, by = "rep.coord")
selected_MT_fpkm <- 
  selected_MT_fpkm %>%
  mutate(repName = gsub("_Mm", "", x = repName)) %>%
  filter(!grepl("MTE", repName)) %>%
  arrange(repName, desc(s_GV.WE))

# write data.frame as .csv
write.csv(x = selected_MT_fpkm, file = "MT_ordered_GV_fpkm.csv")

# read data from csv
selected_MT_fpkm <- read.csv("MT_ordered_GV_fpkm.csv", stringsAsFactors = F, row.names = 1)
  
# write .FASTA for each class separately 
sapply(X = unique(selected_MT_fpkm$repName), 
       FUN = function(X) 
         getSeqWriteFASTAquantiles(ranges_MT = selected_MT_fpkm[selected_MT_fpkm$repName == X, ],
                                   file_name = X))

# all MT families together
all_MT_fpkm <- 
  selected_MT_fpkm %>%
  arrange(desc(s_GV.WE))

all_MT_fpkm <- all_MT_fpkm %>%
  mutate(quantile = ntile(s_GV.WE, 3)) %>%
  mutate(quantile = factor(quantile, labels = c("bottom_33", "mid_33", "top_33")))

invisible(sapply(X = unique(all_MT_fpkm$quantile), 
                 FUN = function(X){
                   
                   # filter by expression quantile
                   all_MT_fpkm <- all_MT_fpkm %>%
                     filter(quantile == X)
                   
                   # get sequence
                   ranges_MT_seq <- getSeq(x = Mmusculus,
                                           names = all_MT_fpkm$chr, 
                                           start = as.integer(all_MT_fpkm$start),  
                                           end = as.integer(all_MT_fpkm$end), 
                                           strand = all_MT_fpkm$rep.strand) 
                   
                   # write.fasta
                   ranges_MT_seq <- as.character(ranges_MT_seq)
                   write.fasta(as.list(ranges_MT_seq),
                               nbchar = 80, 
                               names = paste0(all_MT_fpkm$chr, ":", 
                                              all_MT_fpkm$start, "-", 
                                              all_MT_fpkm$end, ":",
                                              all_MT_fpkm$rep.strand, "|", 
                                              all_MT_fpkm$repName, "| FPKM = ", 
                                              all_MT_fpkm$s_GV.WE),  
                               as.string = TRUE,
                               file.out = paste0("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/motif_discovery/MT/sequences/all_together", 
                                                 "/MT_all_", X, "_fpkm_GV.fasta"),
                               open = "w")
                 }))        