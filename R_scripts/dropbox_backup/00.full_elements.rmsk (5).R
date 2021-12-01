### INFO: 
### DATE: Thu May 02 16:22:37 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2")

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
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# clean repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.MesAur1.0.20160612.RefSeq.clean.fa.out.gz")

######################################################## READ DATA
# read raw repeatMasker
rmsk_df <- readr::read_delim(file = rmsk_path, delim = "\t")

######################################################## MAIN CODE
### clean repeatMasker
rmsk_tidy <- 
  rmsk_df %>%
  dplyr::filter(repFamily != "Low_complexity",
                repClass != "Simple_repeat") %>%
  dplyr::mutate(id = str_c(repName, "_", seqnames, ":", start, "-", end)) %>% 
  GRanges(.)

# get composition of each rmsk ID
rmsk_id <- 
  rmsk_tidy %>% 
  as.data.table(.) %>% 
  .[, list(repName = str_c(repName, collapse = "/"), 
           repClass = str_c(unique(repClass), collapse = "/"), 
           repFamily = str_c(unique(repFamily), collapse = "/"), 
           strand = str_c(unique(strand), collapse = "/")), by = "rmsk_ID"] %>% 
  .[, rmsk_id := as.character(rmsk_id)] %>% 
  .[]

# split repeats by rmsk ID and get range
rmsk_reduced <-
  split(rmsk_tidy, rmsk_tidy$rmsk_id) %>%
  range(., ignore.strand = T)


### find elements which are overlaped by another element
# find overlaps between whole ranges and individual elements
rmsk_overlaps <- findOverlaps(rmsk_reduced, ignore.strand = T)
rmsk_overlaps <- rmsk_overlaps[!isSelfHit(rmsk_overlaps)]

# get names of elements overlaped by some other element
rmsk_overlaping_names <- 
  names(rmsk_reduced[queryHits(rmsk_overlaps)]) %>% 
  unique(.)

# get names of uninterupted elements (= elements not overlaping any elements)
rmsk_whole_names <- names(rmsk_reduced)[!(names(rmsk_reduced) %in% rmsk_overlaping_names)]


### find elements which are completely within other element
# get elements overlaping other element
rmsk_overlaping <- rmsk_reduced[rmsk_overlaping_names]

# find self overlaps, remove self hits
rmsk_within_overlaps <- findOverlaps(rmsk_overlaping, type = "within", ignore.strand = T)
rmsk_within_overlaps <- rmsk_within_overlaps[!isSelfHit(rmsk_within_overlaps)]

# get elements which are inside other elements
rmsk_within <- rmsk_overlaping[unique(queryHits(rmsk_within_overlaps))]

# get names of elements completely within other element
rmsk_within_names <- 
  rmsk_overlaping[unique(queryHits(rmsk_within_overlaps))] %>% 
  names(.)


# check if some of the "within" elements have another element within them (i. e. if they are interupted)
rmsk_within_overlaps_2 <- findOverlaps(rmsk_within, type = "within", ignore.strand = T)
rmsk_within_overlaps_2 <- rmsk_within_overlaps_2[!isSelfHit(rmsk_within_overlaps_2)]

# get names of within elements which have another within element inside them
rmsk_within_interupted_names <- 
  rmsk_within[unique(subjectHits(rmsk_within_overlaps_2))] %>% 
  names(.)

# get names of elements which are completely within another element and are also not interupted
rmsk_within_names <- setdiff(rmsk_within_names, rmsk_within_interupted_names)
  

### find interupted elements
# get difference between all interupted and elements completely within other element
rmsk_interupted_names <- setdiff(rmsk_overlaping_names, rmsk_within_names)


### extract all elements
# whole elements
rmsk_whole <- 
  rmsk_reduced[names(rmsk_reduced) %in% rmsk_whole_names] %>% 
  unlist(.)
rmsk_whole$rmsk_id <- names(rmsk_whole)
rmsk_whole %<>% 
  as.data.frame(.) %>% 
  as.data.table(.) %>% 
  .[rmsk_id, `:=`(strand = i.strand, 
                  repName = i.repName,
                  repClass = i.repClass, 
                  repFamily = i.repFamily), 
    on = "rmsk_id"] %>% 
  .[, insertion_class := "whole"] %>% 
  .[]

# within other elements
rmsk_within <- 
  rmsk_reduced[names(rmsk_reduced) %in% rmsk_within_names] %>% 
  unlist(.)
rmsk_within$rmsk_id <- names(rmsk_within)
rmsk_within %<>% 
  as.data.frame(.) %>% 
  as.data.table(.) %>% 
  .[rmsk_id, `:=`(strand = i.strand, 
                  repName = i.repName,
                  repClass = i.repClass, 
                  repFamily = i.repFamily), 
    on = "rmsk_id"] %>% 
  .[, insertion_class := "within"] %>% 
  .[]

# interupted elements
rmks_interupted <- 
  rmsk_reduced[names(rmsk_reduced) %in% rmsk_interupted_names] %>% 
  unlist(.)
rmks_interupted$rmsk_id <- names(rmks_interupted)
rmks_interupted %<>% 
  as.data.frame(.) %>% 
  as.data.table(.) %>% 
  .[rmsk_id, `:=`(strand = i.strand, 
                  repName = i.repName,
                  repClass = i.repClass, 
                  repFamily = i.repFamily), 
    on = "rmsk_id"] %>% 
  .[, insertion_class := "interupted"] %>% 
  .[]

# join all elements, clean and save
rmsk_all <- 
  rbind(rmsk_whole, rmsk_within, rmks_interupted) %>%
  as_tibble(.) %>% 
  dplyr::select(-width) %>% 
  dplyr::select(seqnames:strand, repName:insertion_class, rmsk_id) %>% 
  dplyr::mutate(rmsk_id = as.integer(rmsk_id)) %>% 
  dplyr::arrange(rmsk_id) %T>%
  readr::write_delim(., path = file.path(genome_dir, basename(rmsk_path) %>% str_replace(., "clean.fa.out.gz", "joined_rmsk_id.fa.out.gz")), delim = "\t")

