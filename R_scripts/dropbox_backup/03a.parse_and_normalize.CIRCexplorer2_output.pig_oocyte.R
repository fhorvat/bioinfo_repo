### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/circRNA_detection/datasets/pig.susScr11/CIRCexplorer2")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# counts table path
counts_tb_path <- list.files(inpath, ".*\\.circularRNA_known\\.txt$", full.names = T)

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/Mapped/pig.susScr11"

# find stats and tracks table
sample_tb_path <- list.files(mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

######################################################## READ DATA
# read counts table
counts_tb <- purrr::map(counts_tb_path, function(path){
  
  # read CIRCExplorer2 output
  circ_tb <- 
    readr::read_delim(path, delim = "\t", 
                      col_names = c("seqnames", "start", "end", "name", "strand", 
                                    "score", "thickStart", "thickEnd", "itemRgb",
                                    "exonCount", "exonSizes", "exonOffset", "readNumber",
                                    "circType", "geneName", "isoformName", "index", "flankIntron")) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove(., "\\.circularRNA_known\\.txt$"))
  
  # return
  return(circ_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

######################################################## MAIN CODE
# clean sample table
sample_tb_clean <- 
  sample_tb %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.PE"), 
                library_size = (library_size / 10e5))

# normalize count table, calculate mean value
count_tb_norm <- 
  counts_tb %>% 
  dplyr::select(sample_id, seqnames, start, end, geneName, count = readNumber) %>% 
  tidyr::unite(col = "coordinates", seqnames, start, end, sep = " ") %>% 
  dplyr::left_join(., sample_tb_clean, by = "sample_id") %>% 
  dplyr::mutate(CPM = round((count / library_size), 3)) %>% 
  dplyr::group_by(coordinates, geneName) %>% 
  dplyr::summarise(CPM = mean(CPM), 
                   read_count = mean(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(desc(CPM)) 

# save
readr::write_csv(count_tb_norm, file.path(outpath, "pig_oocyte.circRNA.CIRCexplorer2.CPM.csv"))



