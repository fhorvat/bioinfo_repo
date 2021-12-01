### INFO: 
### DATE: 08. 09. 2017.  
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_September2017/smallRNA_classification_and_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tibble)

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)

######################################################## PATH VARIABLES
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Split_September2017/smallRNA_classification_and_expression"

bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Tam_2008_Nature_GSE10364/Data/Mapped/STAR_mm10"
bam_list <- list.files(path = bam_path, pattern = "*bam$", full.names = T)
log_list <- list.files(path = bam_path, pattern = "*Log.final.out$", full.names = T)

rmsk_path <- "/common/WORK/fhorvat/reference/mouse/mm10/UCSC/UCSC_repeatMasker_withLTRsubclass.txt.gz"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "GffToGRanges.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
# sample table
# sample_table <- 
#   tibble(sample = str_replace_all(bam_list, "\\/.*\\/|_Aligned.sortedByCoord.out.bam", ""), 
#          log_path = log_list) %>% 
#   dplyr::rowwise() %>%
#   dplyr::mutate(library_size = readr::read_lines(file = log_path) %>% 
#                   magrittr::extract(9) %>% 
#                   str_extract(., "[0-9].*") %>% 
#                   as.integer(.)) %>% 
#   dplyr::select(-log_path)

# repeatMasker
rmsk_gr <-
  read_delim(file = rmsk_path, delim = "\t") %>%
  GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T)

# repeatMasker - rRNA
rmsk_rRNA <- rmsk_gr[rmsk_gr$repClass == "rRNA"]

######################################################## MAIN CODE
### read coordinates of mapped reads from .bam
# create list of genomic coordinates 21-23 bp long

####
bam_file <- bam_list[str_detect(bam_list, "s_GV_19to30")]

# chr4 117159633 117182624
# chr15 81704248 81729919
which_ranges <- GRanges("chr4", IRanges(117159633, 117182624))
which_ranges <- GRanges("chr15", IRanges(81704248, 81729919))
####

# read bam file as data frame
param <- Rsamtools::ScanBamParam(what = c("qname", "rname", "strand", "pos", "seq", "cigar"), 
                                 which = which_ranges, 
                                 tagFilter = list(nM = 0))
bam_in <- Rsamtools::scanBam(bam_file, param = param)

# get aligned sequence
seq_aligned_width <-
  GenomicAlignments::sequenceLayer(x = bam_in[[1]]$seq, cigar = bam_in[[1]]$cigar) %>%
  width(.)

# convert to data.frame, filter for miRNA clusters
bam_in[[1]]$seq <- NULL
bam_in[[1]]$cigar <- NULL
bam_gr <-
  dplyr::bind_cols(bam_in) %>%
  dplyr::mutate(seq_width = seq_aligned_width) %>%
  dplyr::mutate(end = (pos + seq_width - 1),
                full_pos = str_c(rname, ":", pos, "-", end, ":", strand)) %>%
  dplyr::select(-seq_width, seqnames = rname, start = pos, end, strand, full_pos) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  .[width(.) >= 21 & width(.) <= 23]
