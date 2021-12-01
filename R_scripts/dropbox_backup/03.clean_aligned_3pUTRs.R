### INFO: 
### DATE: Sun Nov 04 19:08:17 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/multiz60way.3primeUTRs")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# multiple alignment of 3' UTRs between mouse and rat/dog/human
aligned_UTRs_path <- list.files(inpath, ".*3primeUTRs.txt", full.names = T)

######################################################## READ DATA
# read aligned UTRs
aligned_UTRs <- purrr::map(aligned_UTRs_path, function(path){
  
  # get target name
  target_animal <- str_extract(path, "rn5|canFam3|hg19")
  
  # fread table
  aligned_dt <- fread(path, col.names = c("mouse_coords_seq", str_c(target_animal, "_coords"), str_c(target_animal, "_seq")))
  setkey(aligned_dt, "mouse_coords_seq")
  
}) %>% 
  purrr::reduce(., merge, all = T) %>% 
  .[complete.cases(.)]

######################################################## MAIN CODE
# split mouse coordinates and sequences to 2 columns
aligned_UTRs[, `:=`(mouse_coords = str_remove(mouse_coords_seq, "=.*"), 
                    mouse_seq = str_remove(mouse_coords_seq, ".*="), 
                    mouse_coords_seq = NULL)]

# set column order
setcolorder(aligned_UTRs, c("mouse_coords", "mouse_seq", "rn5_coords", "hg19_coords", "canFam3_coords", "rn5_seq", "hg19_seq", "canFam3_seq"))

#####
## create DNAStringSetsList
aligned_3pUTRs_list <- DNAStringSetList(mm10 = DNAStringSet(x = set_names(aligned_UTRs$mouse_seq, aligned_UTRs$mouse_coords)), 
                                        rn5 = DNAStringSet(x = set_names(aligned_UTRs$rn5_seq, aligned_UTRs$mouse_coords)), 
                                        hg19 = DNAStringSet(x = set_names(aligned_UTRs$hg19_seq, aligned_UTRs$mouse_coords)), 
                                        canFam3 = DNAStringSet(x = set_names(aligned_UTRs$canFam3_seq, aligned_UTRs$mouse_coords)))

# save RDS
saveRDS(aligned_3pUTRs_list, file = file.path(outpath, "aligned_3pUTRs.mm10.rn5_hg19_canFam3.RDS"))

