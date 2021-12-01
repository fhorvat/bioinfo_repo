### INFO: 
### DATE: Mon Oct 28 18:43:00 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/homework")

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
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

#  rmsk coordinates path
rmks_path <- file.path(inpath, "rmsk.mm10.20180919.chr11.fa.out.csv")

######################################################## READ DATA
# read LTRs 
rmsk_tb <- readr::read_csv("rmsk.mm10.20180919.chr11.fa.out.csv")

######################################################## MAIN CODE
### a)
# get LINE1s from RepeatMasker table
line1_tb <- 
  rmsk_tb %>% 
  dplyr::filter(repFamily  == "L1")

# get LINE1s as GRanges
# line with elementNROWS() will remove insertions with fragments on opposite strands
line1_gr <- 
  line1_tb %>% 
  GRanges(.) %>% 
  split(., .$rmsk_id) %>% 
  range(.) %>% 
  .[elementNROWS(.) == 1] %>% 
  unlist(.)
mcols(line1_gr)$rmsk_id <- names(line1_gr)
names(line1_gr) <- NULL

# # get repNames for each LINE1 element
# line1_repnames <- 
#   line1_tb %>% 
#   dplyr::group_by(rmsk_id) %>% 
#   dplyr::summarise(repName = str_c(repName, collapse = "/"), 
#                    strand = str_c(unique(strand), collapse = "/")) %>% 
#   dplyr::mutate(strand = ifelse(strand %in% c("+", "-"), strand, "*"))


### b)
# get ranges of all repeats, remove Low_complexity and Simple_repeat
rmsk_gr <- 
  rmsk_tb %>% 
  dplyr::filter(!(repClass %in% c("Low_complexity", "Simple_repeat"))) %>% 
  GRanges(.)

# find overlaps between LINE1s and other repeats
overlaps <- findOverlaps(line1_gr, rmsk_gr, ignore.strand = T)

# get ID's of LINE1s and overlapping repeats, keep only rmsk_id for ranges which don't overlap any ranges beside their own
overlap_tb <- 
  tibble(rmsk_id = mcols(line1_gr[queryHits(overlaps)])$rmsk_id, 
         insertion_rmsk_id = mcols(rmsk_gr[subjectHits(overlaps)])$rmsk_id) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(insertion_rmsk_id = str_c(unique(insertion_rmsk_id), collapse = "|")) %>%
  dplyr::filter(rmsk_id == insertion_rmsk_id)

# keep only LINE1s without any interuptions
line1_without_interuption <- line1_gr[mcols(line1_gr)$rmsk_id %in% overlap_tb$rmsk_id]


### c)
# find overlaps between all LINE1 insertions and uninterupted
overlaps <- findOverlaps(line1_without_interuption, line1_gr, ignore.strand = T, type = "within")

# get rmsk ID's to table
overlap_tb <- 
  tibble(rmsk_id = mcols(line1_without_interuption[queryHits(overlaps)])$rmsk_id, 
         insertion_rmsk_id = mcols(line1_gr[subjectHits(overlaps)])$rmsk_id) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(insertion_rmsk_id = str_c(unique(insertion_rmsk_id), collapse = "|")) %>%
  dplyr::filter(rmsk_id != insertion_rmsk_id)

# add info about insertion to metadata
mcols(line1_without_interuption)$insert_type <- ifelse(test = mcols(line1_without_interuption)$rmsk_id %in% overlap_tb$rmsk_id, 
                                                       yes = "within", 
                                                       no = "solo")

### d)
# get only full length elements
line1_full_length <- line1_without_interuption[width(line1_without_interuption) >= 5500 & width(line1_without_interuption) <= 6500]

### e)
# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Mmusculus.UCSC.mm10, names = line1_full_length)
names(line1_seq) <- mcols(line1_full_length)$rmsk_id

# find ORFs
line1_orfs <- predORF(line1_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)

# get tidy table
line1_orfs_tb<- 
  line1_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(rmsk_id = seqnames, width) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(longest_orf_1 >= 3600, 
                longest_orf_2 >= 1100)

# filter LINE1 list
line1_final <- line1_full_length[mcols(line1_full_length)$rmsk_id %in% line1_orfs_tb$rmsk_id]
