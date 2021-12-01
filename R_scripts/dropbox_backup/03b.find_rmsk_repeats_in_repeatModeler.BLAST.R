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

# chosen LTR classes path
classes_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate"
classes_tb_path <- file.path(classes_tb_path, "all_LTR_classes 200730.xlsx")

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"


### repeatMasker
# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

# DFAM annotation path
dfam_path <- "/common/WORK/fhorvat/programi/RepeatMasker/RepeatMasker-open-4-0-9-p2/Libraries"
dfam_path <- file.path(dfam_path, "Dfam.embl")

### repeatModeler
# joined repeatModeler path
rmod_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "rmsk.Siomi.20200728.RepeatModeler.joined_rmsk_id.fa.out.gz")

# clean repeatModeler path
rmod_clean_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "rmsk.Siomi.20200728.RepeatModeler.clean.fa.out.gz")


### blast hits path
blast_path <- file.path(inpath, "hamster.sequel.draft-20200302.arrow-families.blastn.txt")

######################################################## READ DATA
# # read joined repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# # read joined repeatModeler
# rmod_tb <- readr::read_delim(rmod_path, delim = "\t")

# read clean repeatModeler
rmod_clean <- readr::read_delim(rmod_clean_path, delim = "\t")

# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path) %>% as_tibble(.)

# read blast hits
blast_results <- readr::read_delim(file = blast_path,
                                   delim = "\t",
                                   col_names = c("query_id", "subject_id", "identity_perc", "alignment_length",
                                                 "mismatches", "gap_open",
                                                 "query_start", "query_end", "subject_start", "subject_end",
                                                 "e_value", "bit_score"))

# read DFAM annotation
dfam_tb <- readr::read_lines(dfam_path)

######################################################## MAIN CODE
### clean data
# clean DFAM annotation
dfam_tb_clean <- 
  tibble(dfam_id = dfam_tb) %>% 
  dplyr::filter(str_detect(dfam_id, "^ID|Type: |SubType: |Species: |^NM")) %>% 
  dplyr::mutate(row_id = cumsum(str_detect(dfam_id, "^ID")), 
                dfam_id = 
                  dfam_id %>% 
                  str_remove_all(., ";.*|CC +") %>% 
                  str_replace(., "^ID +", "ID: ") %>% 
                  str_replace(., "^NM +", "NM: ")) %>% 
  tidyr::separate(dfam_id, into = c("type", "value"), sep = ": ") %>% 
  dplyr::mutate(type = str_replace_all(type, pattern = c("ID" = "dfam_id", 
                                                         "NM" = "repName",
                                                         "^Type$" = "repClass", 
                                                         "SubType" = "repFamily"))) %>% 
  tidyr::pivot_wider(id_cols = row_id, names_from = "type", values_from = "value") %>%
  dplyr::select(dfam_id, repName)

# clean blast results
blast_tb <- 
  blast_results %>% 
  dplyr::mutate(subject_id = str_remove(subject_id, "#.*"), 
                query_id = str_remove(query_id, "#.*")) %>% 
  dplyr::rename(repName_rmod = query_id, dfam_id = subject_id) %>% 
  dplyr::left_join(., dfam_tb_clean, by = "dfam_id") %>% 
  dplyr::select(repName_rmod, dfam_id, repName, everything()) %>% 
  dplyr::mutate(repName = ifelse(is.na(repName), dfam_id, repName))%>% 
  dplyr::group_by(repName_rmod) %>% 
  dplyr::slice_max(alignment_length, n = 1) %>% 
  dplyr::slice_max(identity_perc, n = 1) %>% 
  dplyr::slice_max(bit_score, n = 1, with_ties = F) %>% 
  dplyr::ungroup(.)


### summarize RepeatMasker
# get chosen repeats
classes_tb_filt <- 
  classes_tb %>% 
  dplyr::filter(category_I %in% c("MYSERVx", "MuLV", "IAP", "MMERVKx", "RLTR3x ERV1", "MT", "ORR1")) %>% 
  dplyr::select(repName, category_I)
  
# summarize per category
results_sum <- 
  blast_tb %>% 
  dplyr::right_join(., classes_tb_filt, by = "repName") %>% 
  dplyr::filter(!is.na(repName_rmod)) %>% 
  dplyr::group_by(category_I) %>% 
  dplyr::summarise(repName_rmod = str_c(repName_rmod, collapse = ", "))

