### INFO:
### DATE: Sat Jul 20 21:51:54 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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
genome_dir <- "/common/DB/genome_reference/Muridae/Acomys_cahirinus/AcoCah.GCA_004027535.1"

# clean repeatMasker path
rmsk_path <- list.files(genome_dir, "rmsk\\..*\\.clean\\.fa\\.out\\.gz", full.names = T)

######################################################## READ DATA
# read raw repeatMasker
rmsk_tb <- readr::read_delim(file = rmsk_path, delim = "\t")

######################################################## MAIN CODE
### tidy repeatMasker, convert to GRanges and save as bed
# format as bed, split to repeatClasses
rmsk_bed <-
  rmsk_tb %>%
  dplyr::filter(!(str_detect(repClass, "\\?"))) %>%
  dplyr::mutate(start = start - 1, name = repName, score = 1000, thickStart = start, thickEnd = end, itemRgb = "255,0,0") %>%
  dplyr::select(chrom = seqnames, chromStart = start, chromEnd = end, name, score, strand, thickStart, thickEnd, itemRgb, repClass) %>%
  split(., .$repClass)

# save each repClass as separate .bed
trackDb_all <-
  purrr::map(names(rmsk_bed), function(rep_class){
    
    # subset, save
    rmsk_bed[[rep_class]] %>%
      dplyr::select(-repClass) %>%
      readr::write_delim(., file.path(outpath, rmsk_path %>% basename(.) %>% str_replace(., "\\.clean\\.fa\\.out\\.gz", str_c(".", rep_class, ".bed"))),
                         delim = "\t", col_names = F)
    
    # get line for trackDb.txt
    trackDb <-
      tibble(category = c("track", "parent", "shortLabel", "longLabel",
                          "priority", "spectrum", "maxWindowToDraw", "colorByStrand",
                          "type", "bigDataUrl", ""),
             values = c(str_c("repeatMasker_", rep_class), "repeatMasker", rep_class, str_c(rep_class, " Repeating Elements by RepeatMasker"),
                        "1", "on", "10000000", "50,50,150 150,50,50",
                        "bigBed 6 +",
                        rmsk_path %>% basename(.) %>% str_remove(., "\\.clean\\.fa\\.out\\.gz") %>% str_c(., rep_class, "bb", sep = "."), "")) %>%
      dplyr::mutate(file_out = str_c("\t", category, values, sep = " ")) %$%
      file_out
    
    # return
    return(trackDb)
    
  })

# save repeatMasker trackDb
trackDb_all %>%
  unlist(.) %>%
  readr::write_lines(., path = file.path(outpath, rmsk_path %>% basename(.) %>% str_replace(., "\\.clean\\.fa\\.out\\.gz", ".trackDb")))

