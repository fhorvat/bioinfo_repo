### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/find_hits_in_genome")

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
# library(Biostrings)
# library(BSgenome.Maur.UCSC.Siomi)
# library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bed path
bed_path <- file.path(inpath, "L1_full_length_manual_200707.without_5p_repeat.consensus.BLAST_hits.with_both_ORFs.2kb_upstream.bed")

# RNA-seq path
rnaseq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq/Data/Mapped/STAR_Siomi/5_merged_replicates"
rnaseq_bams_paths <- list.files(rnaseq_path, "\\.bam$", full.names = T)
rnaseq_bams_paths <- rnaseq_bams_paths[1]

# small RNA-seq path
smallrnaseq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Mapped/STAR_Siomi/8_merged_replicates"
smallrnaseq_bams_path <- list.files(rnaseq_path, "\\.bam$", full.names = T)
smallrnaseq_bams_path <- smallrnaseq_bams_path[1]

######################################################## READ DATA
# read bed
bed_gr <- rtracklayer::import.bed(bed_path)

# read data from bam - RNAseq
rnaseq_gr_list <- purrr::map(rnaseq_bams_paths, function(path){
  
  # read bam, 
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = path, use.names = TRUE, param = ScanBamParam(which = bed_gr)) %>% 
    GenomicRanges::grglist(.) %>% 
    GenomicRanges::reduce(., ignore.strand = T) %>% 
    unlist(.)
  
  # return 
  return(bam_gr)
  
}) %>% 
  set_names(., rnaseq_bams_paths %>% basename(.) %>% str_remove(., "\\.bam$"))

######################################################## MAIN CODE
# for each read find to which LINE1 it maps
purrr::map(rnaseq_gr_list, function(gr){
  
  gr <- rnaseq_gr_list[[1]]
  
  # find overlaps
  overlaps <- findOverlaps(gr, bed_gr, ignore.strand = T)
  
  # get all reads and LINE1 hits
  overlaps_tb <- 
    tibble(seqnames = as.character(seqnames(gr[queryHits(overlaps)])), 
           start_read = start(gr[queryHits(overlaps)]), 
           end_read = end(gr[queryHits(overlaps)]), 
           rmsk_id = mcols(bed_gr[subjectHits(overlaps)])$name) %>% 
    tidyr::separate(col = rmsk_id, into = c("hit_coordinates", "rmsk_id"), sep = "\\.") %>% 
    dplyr::mutate(hit_coordinates = str_remove(hit_coordinates, ".*:")) %>% 
    tidyr::separate(col = hit_coordinates, into = c("rmsk_start", "rmsk_end"), sep = "-")
    
  # calculate relative coordinates of hit in repeat
    
})


