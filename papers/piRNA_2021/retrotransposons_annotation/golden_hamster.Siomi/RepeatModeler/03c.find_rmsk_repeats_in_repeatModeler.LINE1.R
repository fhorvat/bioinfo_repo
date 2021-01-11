### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/RepeatModeler")

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


### repeatMasker
# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")


### repeatModeler
# joined repeatModeler path
rmod_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "rmsk.Siomi.20200728.RepeatModeler.joined_rmsk_id.fa.out.gz")

# clean repeatModeler path
rmod_clean_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "rmsk.Siomi.20200728.RepeatModeler.clean.fa.out.gz")


######################################################## READ DATA
# # read joined repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# # read joined repeatModeler
# rmod_tb <- readr::read_delim(rmod_path, delim = "\t")

# read clean repeatModeler
rmod_clean <- readr::read_delim(rmod_clean_path, delim = "\t")

######################################################## MAIN CODE
### clean data
# get RepeatMasker as GRanges
rmsk_gr <- 
  rmsk_clean %>% 
  dplyr::filter(repFamily == "L1") %>%
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), remove = F, sep = " ") %>% 
  GRanges(.)

# get repeatModeler as GRanges
rmod_gr <- 
  rmod_clean %>% 
  dplyr::filter(repFamily == "L1") %>%
  GRanges(.)


### find repeatModeler annotation
# overlap with repeatModeler
overlap <- findOverlaps(rmsk_gr, rmod_gr, ignore.strand = F, minoverlap = 10)

# get repeatMasker and RepeatModeler hits
rmsk_hits <- rmsk_gr[queryHits(overlap)]
rmod_hits <- rmod_gr[subjectHits(overlap)]

# get repeatModeler hits
rmod_class_tb <- 
  tibble(coordinates = rmsk_hits$coordinates, 
         rmsk_repName = rmsk_hits$repName, 
         rmsk_repFamily = rmsk_hits$repFamily, 
         rmsk_width = width(rmsk_hits), 
         rmsk_overlap_width = width(pintersect(rmsk_hits, rmod_hits)), 
         rmod_repName = rmod_hits$repName, 
         rmod_repFamily = rmod_hits$repFamily, 
         rmod_width = width(rmod_hits), 
         rmod_overlap_width = width(pintersect(rmod_hits, rmsk_hits))) %>% 
  dplyr::mutate(rmsk_overlap_percentage = 100*(rmsk_overlap_width / rmsk_width)) %>% 
  dplyr::filter(rmsk_repFamily == rmod_repFamily, 
                rmsk_overlap_percentage > 50) %>%
  dplyr::group_by(rmsk_repFamily, rmod_repName) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::slice_max(count, n = 10) %>% 
  dplyr::group_by(rmsk_repFamily) %>% 
  dplyr::summarise(repeatModeler_repName = str_c(unique(rmod_repName), collapse = ", ")) %T>% 
  readr::write_csv(., file = file.path(outpath, "repeatMasker_repeatModeler.repNames.LINE1.csv"))




