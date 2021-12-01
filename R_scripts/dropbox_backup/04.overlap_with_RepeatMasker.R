### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV/annotation")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

# gaps path
gaps_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hubs/golden_hamster.Siomi/files/gaps/hamster.sequel.draft-20200302.arrow.bed"

# repeatModeler MYSERV annotation path
myserv_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV"
myserv_path <- file.path(myserv_path, "MYSERV.FLI_elements.bed")

######################################################## READ DATA
# # read joined repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# # read gaps
# gaps_gr <- rtracklayer::import(gaps_path)

# read myserv coordinates
myserv_gr <- rtracklayer::import(myserv_path)

######################################################## MAIN CODE
# get LTRs from repeatMasker 
ltr_gr <- 
  rmsk_clean %>% 
  dplyr::filter(repClass == "LTR") %>% 
  GenomicRanges::GRanges(.)

# overlap with repeatMasker 
overlaps <- findOverlaps(myserv_gr, ltr_gr, ignore.strand = F)

# add to table
myserv_tb <- 
  tibble(myserv_name = mcols(myserv_gr)$name[queryHits(overlaps)], 
         repName = mcols(ltr_gr)$repName[subjectHits(overlaps)]) %>% 
  dplyr::group_by(myserv_name) %>% 
  dplyr::summarise(repName = str_c(unique(repName), collapse = ",")) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::count(repName) %>% 
  arrange(-n)

