### INFO: 
### DATE: Thu Aug 22 17:01:47 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

# # set number of threads data.table uses
# nthreads <- 1
# data.table::setDTthreads(nthreads)

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/other/repeats_Maja")

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

library(BSgenome.Mmusculus.UCSC.mm10)
library(stringdist)
library(stringr)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# mkuzman function
getsequences <- function(gr, windowsize = 50, LEN = 150, SOMEVALUE = 55, elementLength = 5000){
  
  # get scanning windows
  seqstart <- start(gr) - elementLength
  chrs <- as.character(seqnames(gr))
  length <- end(gr) + elementLength - seqstart
  starts <- seq(seqstart, seqstart + length, windowsize)
  ends <- starts + LEN
  x <- GRanges(chrs, IRanges(starts, ends))
  
  # extract sequences from genome
  somesequence <- getSeq(BSgenome.Mmusculus.UCSC.mm10, x)
  
  # get distance between all sequences 
  alldists <- as.data.table(stringdistmatrix(somesequence, somesequence))
  alldists[, row := starts]
  alldists <- melt(alldists, id.vars = c("row"))
  alldists[, ":="(col = ends[as.numeric(str_extract(variable, "\\d+"))], 
                  chr = unique(chrs))]
  alldists[value == 0, value := 1000]
  alldists <- alldists[row < col, .(chr, row, col, value)][order(row)]
  
  # get pairs of sequences with string distance lower than cutoff
  x <- alldists[value < SOMEVALUE]
  
  # joins consecutive windows (I guess)?
  x[, follows := 
      {k = (row == dplyr::lag(row) + windowsize);
      l = (col == dplyr::lag(col) + windowsize);
      p = ifelse(k == l & l == 1, 1, 0)
      p = rleid(p)
      ifelse(p %% 2 == 0, p, (p + 1)) %/% 2
      }]
  
  x[, .(chr = unique(chr),
        rowstart = min(row),
        rowend = max(row),
        colchr = unique(chr),
        colstart = min(col),
        colend = max(col),
        value = mean(value)),
    follows]
  
  # x[, .(chr = unique(chr),
  #            start = min(c(row, col)),
  #            end = max(c(row, col)),
  #            value = mean(value)),
  #        follows]
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# full length LINE-1 elements by ZJM 
l1_full_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv"

# annotated uninterupted elements in repeatMasker path
rmsk_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz"

######################################################## READ DATA
# read full length LINE-1 elements
l1_full_tb <- readr::read_csv(l1_full_path)

# read annotated uninterupted elements in repeatMasker path
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
### prepare data
# get MTA whole elements
mta_tb <- 
  rmsk_tb %>% 
  dplyr::filter(insertion_class %in% c("within", "whole")) %>% 
  dplyr::filter(str_detect(repName, "MTA")) %>% 
  dplyr::filter(repName == "MTA_Mm/MTA_Mm-int/MTA_Mm") %T>% 
  readr::write_csv(., file.path(outpath, "MTA_Mm.uninterrupted.csv"))

# get uninterrupted MTA as GRanges
mta_gr <- GRanges(mta_tb)

# apply Maja's function to uninterrupted MTA elements
mta_results <- purrr::map(1:20, function(n){
  
  # get one element and find ranges
  mta_gr[n] %>% 
    getsequences(., windowsize = 50, LEN = 150, SOMEVALUE = 55, elementLength = 1800)
  
})




