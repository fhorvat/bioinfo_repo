### INFO: 
### DATE: Tue Jul 09 16:09:23 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# maf .bed paths
maf_paths <- file.path(inpath, "UCSC_maf")

# blat .bed paths
blat_paths <- file.path(inpath, "blat_results")

# list maf .bed files
maf_bed_paths <- list.files(maf_paths, pattern = "\\.mm10\\.lnc1_3prime\\.top_result\\.bed", full.names = T)

# blat table path
blat_tb_path <- file.path(blat_paths, "mm10.lnc1_3prime.top_results.blat.csv")

# genomes info path
genomes_info_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/genomes", 
                               "genome_info.Muridae_with_mouse_strains.tidy.20190706.csv")

######################################################## READ DATA
# read maf .bed, join to table
maf_tb <- 
  purrr::map(maf_bed_paths, function(path){
    
    # read bed, transform to table, tidy
    rtracklayer::import.bed(path) %>% 
      as_tibble(.) %>% 
      dplyr::mutate(species = basename(path) %>% str_remove(., "\\.mm10\\.lnc1_3prime\\.top_result\\.bed")) %>% 
      tidyr::unite(coords, seqnames, start, end, sep = " ") %>% 
      dplyr::select(coords.maf = coords, width.maf = width, species)
    
  }) %>% 
  dplyr::bind_rows(.) %>% 
  arrange(species)

# read blat table
blat_tb <- 
  readr::read_csv(blat_tb_path) %>% 
  dplyr::mutate(species = str_remove(species, "_v1$")) %>% 
  dplyr::select(coords.blat = coords, width.blat = width, strand, score.blat = score, perc_identity.blat = perc_identity, 
                mm10_width.blat = mm10_width, species)

# read genomes info
genomes_info <- 
  readr::read_csv(genomes_info_path) %>% 
  dplyr::mutate(species = str_remove(genome_name_short, "_v1$"))
  
######################################################## MAIN CODE
# join .bed tables
bed_tidy <- 
  full_join(blat_tb, maf_tb, by = "species") %>% 
  dplyr::select(species, coords.blat, coords.maf, width.blat, width.maf, mm10_width.blat, everything()) %T>%
  readr::write_csv(., file.path(outpath, "results", "mm10.lnc1_3prime.top_results.blat_vs_maf.csv"))

### manually currated coordinates
## blat 
# blat species
blat_species <- c("129S1_SvImJ", "A_J", "AKR_J", "ApoSyl", "BALB_cJ", "C3H_HeJ", "C57BL_6NJ",
                  "CAST_EiJ", "CBA_J", "DBA_2J", "FVB_NJ", "LP_J", "mm10", "MusCar", "MusPah",
                  "MusSpi", "NOD_ShiLtJ", "NZO_HlLtJ", "SPRET_EiJ", "WSB_EiJ", "PWK_PhJ")

# maf species
maf_species <- c("rn6")

# blat coordinates
blat_bed <- 
  bed_tidy %>% 
  dplyr::filter(species %in% blat_species) %>% 
  dplyr::select(coords.blat, strand, species) %>% 
  tidyr::separate(coords.blat, c("seqnames", "start", "end"), sep = " ")

# maf coordinates
maf_bed <- 
  bed_tidy %>% 
  dplyr::filter(species %in% maf_species) %>% 
  dplyr::select(coords.maf, strand, species) %>% 
  tidyr::separate(coords.maf, c("seqnames", "start", "end"), sep = " ")

# add genomes info to selected species coordinates, write table
dplyr::bind_rows(blat_bed, maf_bed) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::left_join(., genomes_info, by = "species") %>% 
  dplyr::select(species, coordinates, strand, organism_scientific, assembly_name, organism_common) %T>%
  readr::write_csv(., "lnc1_3prime.Muridae.genomes_info.csv")
  
# join to one table, split by species
all_bed <- 
  bind_rows(blat_bed, maf_bed) %>% 
  split(., .$species)

# save as .bed
purrr::map(names(all_bed), function(species){
  
  # get one species, save as .bed
  all_bed[[species]] %>% 
    GRanges(.) %T>% 
    rtracklayer::export.bed(., file.path(outpath, "results", str_c(species, ".lnc1_3prime.currated.bed")))
  
})

