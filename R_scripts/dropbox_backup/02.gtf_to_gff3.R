### INFO:
### DATE: Sat Jul 20 21:51:54 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Zuzka/2019/hamster_LINE1_expression/ensembl_annotation")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# ensembl gtf path
ensembl_path <- file.path(genome_dir, "ensembl.93.MesAur1.0.20180920.gtf.gz")

######################################################## READ DATA
# read ensembl gtf
ensembl_tb <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# get exon coordinates from .gtf
ensembl_gr <-
  ensembl_tb %>%
  gtfToGRanges(., filter = "exon")

# save as gff3
ensembl_genes <- 
  ensembl_tb %>% 
  gtfToGRanges(.) %T>% 
  rtracklayer::export.gff3(., file.path(outpath, ensembl_path %>% basename(.) %>% str_replace(., "\\.gtf\\.gz", ".gff3")))


# # save each repClass as separate .bed
# trackDb_all <-
#   purrr::map(names(rmsk_bed), function(rep_class){
#     
#     # subset, save
#     rmsk_bed[[rep_class]] %>%
#       dplyr::select(-repClass) %>%
#       readr::write_delim(., file.path(outpath, rmsk_path %>% basename(.) %>% str_replace(., "\\.clean\\.fa\\.out\\.gz", str_c(".", rep_class, ".bed"))),
#                          delim = "\t", col_names = F)
#     
#     # get line for trackDb.txt
#     trackDb <-
#       tibble(category = c("track", "parent", "shortLabel", "longLabel",
#                           "priority", "spectrum", "maxWindowToDraw", "colorByStrand",
#                           "type", "bigDataUrl", ""),
#              values = c(str_c("repeatMasker_", rep_class), "repeatMasker", rep_class, str_c(rep_class, " Repeating Elements by RepeatMasker"),
#                         "1", "on", "10000000", "50,50,150 150,50,50",
#                         "bigBed 6 +",
#                         rmsk_path %>% basename(.) %>% str_remove(., "\\.clean\\.fa\\.out\\.gz") %>% str_c(., rep_class, "bb", sep = "."), "")) %>%
#       dplyr::mutate(file_out = str_c("\t", category, values, sep = " ")) %$%
#       file_out
#     
#     # return
#     return(trackDb)
#     
#   })
# 
# # save repeatMasker trackDb
# trackDb_all %>%
#   unlist(.) %>%
#   readr::write_lines(., path = file.path(outpath, rmsk_path %>% basename(.) %>% str_replace(., "\\.clean\\.fa\\.out\\.gz", ".trackDb")))
