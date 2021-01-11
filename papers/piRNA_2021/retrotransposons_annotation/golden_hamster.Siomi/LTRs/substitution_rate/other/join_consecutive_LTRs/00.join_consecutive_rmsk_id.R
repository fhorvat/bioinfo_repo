### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate")

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

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

# raw repeatMasker
rmsk_raw_path <- file.path(genome_dir, "rmsk.Siomi.20200701.raw.fa.out.gz")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# read raw repeatMasker
rmsk_raw <- readr::read_table2(file = rmsk_raw_path, skip = 3, 
                               col_names = c("bit_score", "perc_substitution", "perc_deletion", "perc_insertion", 
                                             "seqnames", "query_start", "query_end", "query_left", "strand",
                                             "repName", "repClass_repFamily", "repeat_begin", "repeat_start", "repeat_end", 
                                             "rmsk_id", "tmp"))

######################################################## MAIN CODE
### filter raw repeatMasker table
# join consecutive inserts with same rmsk_id
rmsk_joined <- 
  rmsk_raw %>% 
  dplyr::select(seqnames, start = query_start, end = query_end, strand, 
                bit_score, perc_substitution, perc_deletion, perc_insertion,
                repName, repClass_repFamily, 
                repeat_begin, repeat_start, repeat_end, 
                rmsk_id) %>% 
  dplyr::mutate(strand = str_replace(strand, "C", "-")) %>% 
  dplyr::mutate(repeatStart = ifelse(strand == "+", repeat_begin, repeat_end), 
                repeatEnd = repeat_start, 
                repeatLeft = ifelse(str_detect(repeat_begin, "\\("), repeat_begin, repeat_end)) %>% 
  dplyr::mutate(repeatStart = as.numeric(repeatStart),
                repeatEnd = as.numeric(repeatEnd), 
                repeatLeft = str_remove_all(repeatLeft, "\\(|\\)") %>% as.numeric(.)) %>% 
  dplyr::select(-c(repeat_begin, repeat_start, repeat_end)) %>% 
  dplyr::select(seqnames:repClass_repFamily, repeatStart:repeatLeft, rmsk_id) 
                
# continue
rmsk_joined_2 <- 
  rmsk_joined %>% 
  dplyr::mutate(repeatWidth = repeatEnd - repeatStart + 1) %>% 
  dplyr::mutate(count_substitution = round((perc_substitution / 100) * repeatWidth), 
                count_deletion = round((perc_deletion / 100) * repeatWidth), 
                count_insertion = round((perc_insertion / 100) * repeatWidth)) %>% 
  dplyr::select(-c(perc_substitution, perc_deletion, perc_insertion)) %>% 
  dplyr::mutate(rmsk_id_consecutive = str_c(rmsk_id, data.table::rleid(rmsk_id), sep = ".")) %>%
  dplyr::group_by(rmsk_id_consecutive) %>%
  dplyr::summarize(seqnames = unique(seqnames),
                   start = min(start),
                   end = max(end),
                   strand = str_c(unique(strand), collapse = "/"),
                   bit_score = sum(bit_score), 
                   repName = str_c(repName, collapse = "/"),
                   repClass_repFamily = str_c(unique(repClass_repFamily), collapse = "|"),
                   repeatStart = min(repeatStart), 
                   repeatEnd = max(repeatEnd), 
                   repeatLeft = min(repeatLeft),
                   repeatWidth = sum(repeatWidth), 
                   count_substitution = sum(count_substitution),
                   count_deletion = sum(count_deletion), 
                   count_insertion = sum(count_deletion),
                   rmsk_id = unique(rmsk_id)) %>%
  dplyr::mutate(perc_substitution = 100 * round(count_substitution / repeatWidth, 3),
                perc_deletion = 100 * round(count_deletion / repeatWidth, 3), 
                perc_insertion = 100 * round(count_insertion / repeatWidth, 3)) %>% 
  dplyr::select(seqnames:bit_score, repName, repClass_repFamily, repeatStart:repeatLeft, repeatWidth, perc_substitution:perc_insertion, rmsk_id) %>%
  arrange(start)

# continue
rmsk_joined_3 <- 
  rmsk_joined_2 %>% 
  dplyr::mutate(repeatAbsWidth = repeatEnd - repeatStart + 1) %>% 
  dplyr::select(seqnames:bit_score, repName, repClass_repFamily, repeatStart:repeatLeft, repeatWidth, repeatAbsWidth, perc_substitution:perc_insertion, rmsk_id) %>% 
  dplyr::arrange(seqnames, start)

# save
readr::write_delim(rmsk_joined_3, file.path(outpath, "rmsk.Siomi.20200701.joined_consecutive.fa.out.gz"), delim = "\t")
