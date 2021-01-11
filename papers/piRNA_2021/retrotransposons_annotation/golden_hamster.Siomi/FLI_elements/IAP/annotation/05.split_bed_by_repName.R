### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/annotation")

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

library(BSgenome.Maur.UCSC.Siomi)
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(stringdist)
library(msa)
library(DECIPHER)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bed path
bed_path <- file.path(inpath, "IAP.potentially_young.ordered_by_ORFs.20201031.top_110.bed")

######################################################## READ DATA
# read bed
bed_gr <- rtracklayer::import.bed(bed_path)

######################################################## MAIN CODE
# split by repName
mcols(bed_gr)$repName <- str_remove(mcols(bed_gr)$name, "\\..*")
names(bed_gr) <- mcols(bed_gr)$name
bed_gr_list <- split(bed_gr, mcols(bed_gr)$repName)

# save separately
purrr::map(names(bed_gr_list), function(rmsk_name){
  
  bed_gr_list[[rmsk_name]] %>% 
    rtracklayer::export.bed(., file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031", rmsk_name, "bed", sep = ".")))
  
})
