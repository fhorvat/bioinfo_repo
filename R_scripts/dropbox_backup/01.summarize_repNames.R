### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate")

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

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# repeatMasker annotation path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# read repeatMasker table
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# count LTRs by repName
ltrs_tb <- 
  rmsk_tb %>% 
  dplyr::filter(repClass == "LTR") %>% 
  dplyr::group_by(repName) %>% 
  dplyr::summarise(repFamily = unique(repFamily), 
                   count = n())

# save
readr::write_csv(ltrs_tb, file.path(outpath, "all_LTR_classes.csv"))


### get associations
# get LTRs
ltrs_tb <- 
  rmsk_tb %>% 
  dplyr::filter(repClass == "LTR")

# LTR name
ltr_name <- "MuLV-int"

# filter
ltr_assotiation <- 
  ltrs_tb %>% 
  dplyr::filter(str_detect(repName, ltr_name)) %$%
  repName %>% 
  str_remove_all(., ltr_name) %>% 
  str_replace_all(., "/", " ") %>% 
  str_trim(.) %>% 
  .[. != ""] %>% 
  str_split(., " ") %>% 
  unlist(.) %>% 
  .[. != ""] %>% 
  tibble(repName = .) %>% 
  dplyr::count(repName) %>% 
  dplyr::arrange(desc(n))
  
