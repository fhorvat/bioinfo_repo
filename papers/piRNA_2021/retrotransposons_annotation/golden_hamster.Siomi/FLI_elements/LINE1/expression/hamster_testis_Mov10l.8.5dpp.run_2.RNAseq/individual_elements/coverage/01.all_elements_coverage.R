### INFO: finds coverage of repetitive elements from list and outputs total coverage ratio and longest coverage stretch in each element
### DATE: Sun Dec 13 20:30:06 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/expression/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq/individual_elements/coverage")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# bigWigs path
bw_path <- "../../bw_subset"

# bed path
bed_path <- "../../../.."
bed_path <- list.files(bed_path, ".*\\.FLI_elements\\.bed", full.names = T)

# list bigWig files
bw_path_list <- list.files(bw_path, pattern = "\\.bw$", full.names = T)
bw_path_list <- bw_path_list[!str_detect(bw_path_list, "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE.bw")]

######################################################## READ DATA
# read .bed
rmsk_bed <- rtracklayer::import.bed(bed_path)

# read coverage from bigWig
element_coverage_list <- purrr::map(bw_path_list, function(path){
  
  # read bigWig file
  coverage <- rtracklayer::import(path)
  
  # add genotype 
  mcols(coverage)$genotype <- str_extract(basename(path), "Mov10l1_WT|Mov10l1_KO")
  
  # return
  return(coverage)
  
})


######################################################## MAIN CODE
### prepare files
# create element table
rmsk_tb <- 
  rmsk_bed %>% 
  as_tibble(.)

# join coverage for each genotype
element_coverage_list <- do.call(c, element_coverage_list)
element_coverage_list <- split(element_coverage_list, mcols(element_coverage_list)$genotype)


### find coverage for each element
element_coverage_tb <- purrr::map(names(element_coverage_list), function(genotype){
  
  # get one element
  element_coverage_full <- element_coverage_list[[genotype]]
  
  # reduce coverage
  element_coverage_full <- reduce(element_coverage_full, ignore.strand = T)
  
  # find overlaps between two GRanges
  hits <- findOverlaps(rmsk_bed, element_coverage_full, ignore.strand = T)
  
  # extract all overalaping features from subject as list
  element_coverage <- extractList(element_coverage_full, as(hits, "List"))
  
  # intersect with element coordinates
  element_coverage <- pintersect(element_coverage, rmsk_bed, ignore.strand = T)
  
  # set names
  names(element_coverage) <- mcols(rmsk_bed)$name 
  
  
  ### find longest stretch and percentage of coverage for each element
  # unlist
  element_coverage_tb <- 
    element_coverage %>% 
    unlist(.)
  mcols(element_coverage_tb)$name <- names(element_coverage_tb)
  
  # convert to table, summarize 
  element_coverage_sum <- 
    element_coverage_tb %>% 
    as_tibble(.) %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise(coverage_total = sum(width), 
                     coverage_longest = max(width)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::left_join(., rmsk_tb %>% dplyr::select(name, width), by = "name") %>% 
    dplyr::mutate(total_coverage_ratio = round((coverage_total / width), 3), 
                  longest_coverage_ratio = round((coverage_longest / width), 3)) %>% 
    dplyr::select(name, total_coverage_ratio, longest_coverage_ratio) %>% 
    tidyr::pivot_longer(-name, names_to = "ratio_name", values_to = "ratio_value") %>% 
    dplyr::mutate(genotype = genotype)

  # return 
  return(element_coverage_sum)
  
}) %>% 
  dplyr::bind_rows(.)

# join to one table
element_coverage_wide <-
  element_coverage_tb %>% 
  dplyr::filter(ratio_name == "longest_coverage_ratio") %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_KO"))) %>% 
  dplyr::arrange(genotype) %>% 
  tidyr::pivot_wider(id_cols = name, names_from = c("ratio_name", "genotype"), values_from = "ratio_value", names_sep = ".")

# save
readr::write_csv(element_coverage_wide, file.path(outpath, str_replace(basename(bed_path), "\\.bed", ".coverage_ratio.csv")))
