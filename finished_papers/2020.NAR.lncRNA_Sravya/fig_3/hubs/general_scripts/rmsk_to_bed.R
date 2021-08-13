### INFO: 
### DATE: Sat Jul 20 21:51:54 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_locus/results/repeatMasker")

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

# list all repeatMasker out files
rmsk_files <- 
  list.files(inpath, pattern = "*.\\.fa\\.out", recursive = T, full.names = T) %>% 
  .[!str_detect(., "backup")]

######################################################## READ DATA

######################################################## MAIN CODE
# read and tidy repeatMasker, convert to GRanges and save as bed
purrr::map(rmsk_files, function(rmsk_path){
  
  cat("\n", rmsk_path, "\n", "\n")
  
  # read
  rmsk_tb <- 
    readr::read_table2(file = rmsk_path, skip = 3, col_names = F) %>%
    dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_id = X15) %>%
    tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
    dplyr::mutate(strand = replace(strand, strand == "C", "-"))
  
  # format as bed, split to repeatClasses
  rmsk_bed <- 
    rmsk_tb %>%  
    dplyr::mutate(start = start - 1, name = repName, score = 1000, thickStart = start, thickEnd = end, itemRgb = "255,0,0") %>% 
    dplyr::select(chrom = seqnames, chromStart = start, chromEnd = end, name, score, strand, thickStart, thickEnd, itemRgb, repClass) %>% 
    split(., .$repClass)
  
  # save each repClass as separate .bed
  trackDb_all <- 
    purrr::map(names(rmsk_bed), function(rep_class){
    
    # subset, save
    rmsk_bed[[rep_class]] %>%
      dplyr::select(-repClass) %>% 
      readr::write_delim(., str_replace(rmsk_path, "fa.out", str_c(rep_class, ".bed")), delim = "\t", col_names = F)
    
    # write trackDb.txt for UCSC hub for repeatMasker groups
    trackDb <- 
      tibble(category = c("track", "parent", "shortLabel", "longLabel", 
                          "priority", "spectrum", "maxWindowToDraw", "colorByStrand", 
                          "type", "bigDataUrl", ""), 
             values = c(str_c("repeatMasker_", rep_class), "repeatMasker", rep_class, str_c(rep_class, " Repeating Elements by RepeatMasker"), 
                        "1", "on", "10000000", "50,50,150 150,50,50", 
                        "bigBed 6 +", str_c((basename(rmsk_path) %>% str_remove(., "\\.fa\\.out")), rep_class, "bb", sep = "."), "")) %>% 
      dplyr::mutate(file_out = str_c("\t", category, values, sep = " ")) %$% 
      file_out
    
    # return 
    return(trackDb)
    
  })
  
  # save repeatMasker trackDb
  trackDb_all %>% 
    unlist(.) %>% 
    readr::write_lines(., path = str_replace(rmsk_path, "out", "trackDb"))
  
  # # format to gff3, save
  # rmsk_tb %>%  
  #   GRanges(.) %T>% 
  #   rtracklayer::export.gff3(., file.path(str_replace(rmsk_path, "\\.fa\\.out", ".gff3")))
  
})
