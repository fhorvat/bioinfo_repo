### INFO: 
### DATE: Thu May 02 16:22:37 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly")

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

# raw repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.raw.fa.out.gz")

# raw repeatModeler path
rmod_path <- file.path(genome_dir, "RepeatModeler/RepeatMasker", "hamster.fasta.out")

# LINE1 full length hits 
line1_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/expression",
                        "L1_full_length_manual_200707.consensus.ORFs_blast_hits.saf")

######################################################## READ DATA
# read raw repeatMasker
rmsk_df <- readr::read_table2(file = rmsk_path, skip = 3, col_names = F)

# read raw repeatModeler annotation mapped with repeatMasker
rmod_df <- readr::read_table2(file = rmod_path, skip = 3, col_names = F)

# read LINE1 full length insertions
line1_tb <- readr::read_delim(line1_path, delim = "\t")

######################################################## MAIN CODE
# clean repeatMasker
rmsk_clean <- 
  rmsk_df %>%
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_id = X15) %>%
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) 

# clean repeatModeler
rmod_clean <- 
  rmod_df %>%
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_id = X15) %>%
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) 


### overlap annotations
# create GRanges
rmsk_gr <- GRanges(rmsk_clean)
rmod_gr <- GRanges(rmod_clean)

# find overlaps between RepeatMasker and RepeatModeler
overlaps <- findOverlaps(rmod_gr, rmod_gr, ignore.strand = T)

# find guys from RepeatModeler which don't overlap anything in RepeatMasker
rmod_only <- rmod_gr[setdiff(1:length(rmod_gr), unique(queryHits(overlaps)))]

# find guys in RepeatMasker annotation which don't overlap anything in RepeatModeler
rmsk_only <- rmsk_gr[setdiff(1:length(rmsk_gr), unique(subjectHits(overlaps)))]


### check full length insertions
# LINE1 insertions to GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::select(seqnames = Chr, start = Start, end = End, strand = Strand, rmsk_id = GeneID) %>% 
  GRanges(.)

# overlap 
overlap <- findOverlaps(line1_gr, rmod_gr, ignore.strand = F)

# is there any LINE1 which was not annotated by RepeatModeler? 
all(1:length(line1_gr) %in% queryHits(overlap))

# for each insertion get all fragments overlapping it
line1_overlaps <-
  tibble(line1_id = line1_gr[queryHits(overlap)]$rmsk_id, 
         line1_width = width(line1_gr[queryHits(overlap)]), 
         hit_id = rmod_gr[subjectHits(overlap)]$rmsk_id) %>% 
  dplyr::left_join(., rmod_clean, by = c("hit_id" = "rmsk_id")) %>% 
  dplyr::group_by(line1_id) %>% 
  dplyr::summarise(line1_width = unique(line1_width), 
                   hit_id = str_c(unique(hit_id), collapse = "|"), 
                   seqnames = unique(seqnames), 
                   start = min(start), 
                   end = max(end), 
                   strand = unique(strand), 
                   repName = str_c(unique(repName), collapse = "|")) %>% 
  dplyr::mutate(rmod_width = end - start + 1)
  


