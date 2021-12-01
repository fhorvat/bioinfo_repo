### INFO:
### DATE: Thu Jul 30 16:34:41 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/expression/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/individual_expression")

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

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# library size path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.deduplicated.smallRNAseq/Data/Mapped/STAR_Siomi.multimappers"
library_size_path <- file.path(mapped_path, "4_library_size/library_sizes.txt")

# bed path
bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE"
bed_path <- file.path(bed_path, "LINE1.FLI_elements.bed")

# bam path
bam_list <- file.path(inpath, "../bam_subset")
bam_list <- list.files(bam_list, "\\.bam$", full.names = T)

######################################################## READ DATA
# read library size table
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read .bed
rmsk_bed <- rtracklayer::import.bed(bed_path)

######################################################## MAIN CODE
### for each alignment find one random sense and antisense hit
# in a loop for each bam
fli_tb <- purrr::map(bam_list, function(bam_path){
  
  # read alignments from .bam
  bam_alignments <- GenomicAlignments::readGAlignmentsList(bam_path, use.names = T)
  
  # find overlaps with repeatMasker for each alignment
  bam_gr <-
    bam_alignments %>%
    unlist(.)
  
  # find all overlaps
  overlaps <- findOverlaps(bam_gr, rmsk_bed, ignore.strand = T)
  
  # get hits
  bam_hits <- bam_gr[queryHits(overlaps)]
  rmsk_hits <- rmsk_bed[subjectHits(overlaps)]
  
  # summarize expression per repName, read width and sense
  overlaps_tb <-
    tibble(read_name = names(bam_hits),
           read_strand = as.character(strand(bam_hits)),
           read_width = width(bam_hits),
           rmsk_name = mcols(rmsk_hits)$name,
           rmsk_strand = as.character(strand(rmsk_hits))) %>%
    dplyr::mutate(sense = ifelse(read_strand == rmsk_strand, "sense", "antisense")) %>%
    dplyr::group_by(read_name, sense) %>%
    dplyr::sample_n(1) %>%
    dplyr::ungroup(.) %>%
    dplyr::group_by(rmsk_name, sense, read_width) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup(.) %>% 
    dplyr::mutate(sample_id = basename(bam_path) %>% str_remove(., "\\.bam$"))
  
}) %>% 
  dplyr::bind_rows(.)

# summarize for sample id
fli_sum <- 
  fli_tb %>% 
  dplyr::group_by(rmsk_name, sample_id, sense) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(library_size = (library_size / 1e6), 
                rpm = count / library_size) %>% 
  tidyr::pivot_wider(id_cols = c(rmsk_name, sample_id), values_from = rpm, names_from = sense) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "WT|KO")) %>% 
  dplyr::group_by(rmsk_name, genotype) %>% 
  dplyr::summarise(sense = mean(sense), 
                   antisense = mean(antisense)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = rmsk_name, names_from = genotype, values_from = c(sense, antisense), values_fill = 0)
