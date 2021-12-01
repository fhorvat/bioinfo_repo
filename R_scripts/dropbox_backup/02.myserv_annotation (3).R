### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV.repeatMasker/annotation")

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

######################################################## READ DATA
# # read joined repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# # read gaps
# gaps_gr <- rtracklayer::import(gaps_path)

######################################################## MAIN CODE
# get widths
rmsk_clean %<>%
  dplyr::mutate(rmsk_id = as.character(rmsk_id),
                width = end - start + 1)

# filter
rmsk_tb_filt <-
  rmsk_clean %>%
  dplyr::filter(repName %in% c("MYSERV6-int")) %>%
  dplyr::filter(width > 4000) %T>%
  readr::write_csv(., file = file.path(outpath, "MYSERV6-int.longer_than_4kb.csv"))

# save as .bed
rmsk_bed <- 
  rmsk_tb_filt %>% 
  GRanges(.)
names(rmsk_bed) <- as.character(mcols(rmsk_bed)$rmsk_id) 
rtracklayer::export.bed(rmsk_bed, file.path(outpath, "MYSERV6-int.longer_than_4kb.bed"))





#### joined RepeatMasker
# # # filter 
# rmsk_tb_filt <- 
#   rmsk_tb %>% 
#   dplyr::filter(str_detect(repName, "MYSERV6-int")) %>% 
#   dplyr::mutate(width = end - start + 1) %>% 
#   dplyr::filter(width > 4000) %>% 
#   dplyr::filter(insertion_class != "interupted") %>% 
#   dplyr::arrange(desc(width))
# 
# ### remove those which overlap gaps in assembly
# # get rmsk_id
# rmsk_gaps <- 
#   subsetByOverlaps(GRanges(rmsk_tb_filt), gaps_gr, ignore.strand = T) %$% 
#   rmsk_id
# 
# # filter
# rmsk_tb_filt %<>% 
#   dplyr::filter(!(rmsk_id %in% rmsk_gaps))

