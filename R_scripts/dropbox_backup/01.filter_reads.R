### INFO: 
### DATE: Tue Oct 01 22:17:29 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/datasets/Dicer_Mili_KO/Mapped/mm10_masked/1_mapped/test")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get the list of the reads
bam_paths <- file.path(inpath, c("s_GV_DBL_old_r1.PE.bam", "s_GV_DBL_old_r1.PE.FH.bam"))

######################################################## READ DATA

######################################################## MAIN CODE
# read bam, make names of read pairs unique
bam_list <- purrr::map(bam_paths, function(bam_path){
  
  bam <- GenomicAlignments::readGAlignmentsList(file = bam_path, use.names = T, param = ScanBamParam(what = c("qname")))
  names(bam) <- make.unique(names(bam))
  names(bam)[!str_detect(names(bam), "\\.[0-9]+$")] <- str_c(names(bam)[!str_detect(names(bam), "\\.[0-9]+$")], ".0")
  return(bam)
  
}) 

# create table with bam alignments, remove all pairs for each at least one read has a deletion, insertion, splicing or soft clip 
bam_filter <- purrr::map(1:length(bam_list), function(n){
  
  bam_dt <-
    bam_list[[n]] %>%
    unlist(.) %>%
    as_tibble(., rownames = "read_name") %>% 
    as.data.table(.) %>% 
    .[, `:=`(read_name = str_remove(read_name, "(?<=\\.[0-9]{1,1000})\\.[1-9]+"))] %>%
    .[.[, .I[!any(str_detect(cigar, "S|N|I|D"))], by = read_name]$V1] %>% 
    .[, `:=`(mate1 = str_c(start[1], "_", start[2])), by = "read_name"] %>% 
    .[.[, .I[1:2], by = mate1]$V1]
  
})

# get both GRanges
zoe_gr <- bam_filter[[1]] %>% GRanges(.)
fh_gr <- bam_filter[[2]] %>% GRanges(.)

# get area GRanges
test_gr <- GRanges(seqnames = "chr18", ranges = IRanges(start = 4222494, end = 4228835))

# get reads in area for both
zoe_gr_test <- subsetByOverlaps(zoe_gr, test_gr, ignore.strand = T)
fh_gr_test <- subsetByOverlaps(fh_gr, test_gr, ignore.strand = T)

# get reads in Zoe not in other
zoe_gr_test[!(zoe_gr_test$qname %in% fh_gr_test$qname)]

fh_gr_test






# create bigWig track
bam_coverage <-
  bam_dt %>%
  GRanges(.) %>%
  coverage(.) %>%
  GRanges(.)

# save as test bigWig
rtracklayer::export(object = bam_coverage,  con = str_replace(bam_path, "\\.bam$", ".test.bw"))



# create 2 tracks and test, if they are the same continue with removing multimappers which have same 1 read in a pair, but different 2nd 

dplyr::group_by(read_name) %>%
  dplyr::mutate(mate1 = str_c(start[1], "_", end[1])) %>%
  dplyr::ungroup(.) %>%
  dplyr::group_by(mate1) %>%
  dplyr::slice(1:2) %>%
  dplyr::ungroup(.)


test_read <- bam_tb[str_detect(bam_tb$read_name, "NB502156:92:HHGGHBGXC:4:12404:2702:19288"), ]
bam_filt[str_detect(names(bam_filt), "NB502156:92:HHGGHBGXC:4:12404:2702:19288")]
