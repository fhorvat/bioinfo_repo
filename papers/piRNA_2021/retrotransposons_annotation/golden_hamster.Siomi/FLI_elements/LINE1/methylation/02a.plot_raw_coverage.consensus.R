### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/methylation/mapped/Nov_2020/FLI_consensus")

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

# mapped path
mapped_path <- inpath

# bigWigs path
bw_path <- list.files(inpath, ".*raw_coverage.bw$", full.names = T) 

# sequences path
seq_path <- "../../../bismark_index/FLI_consensus"
seq_path <- list.files(seq_path, ".*\\.fasta$", full.names = T)

######################################################## READ DATA
# read data
coverage_tb_full <- purrr::map(bw_path, function(path){
  
  # read bw
  bw_gr <- rtracklayer::import.bw(path, as = "RleList")
  
  # get coverage
  coverage_list <- 
    bw_gr %>% 
    as(., "IntegerList") %>% 
    unlist(.) %>% 
    split(., names(.)) 
  
  # set to table
  coverage_tb <- purrr::map(names(coverage_list), function(seqname){
    
    # create table
    cov_tb <- 
      coverage_list[[seqname]] %>% 
      tibble(coverage = .) %>% 
      dplyr::mutate(repName = seqname,
                    sample_id = basename(path), 
                    pos = 1:nrow(.)) %>% 
      dplyr::select(repName, sample_id, pos, coverage)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# read consensus sequence
cons_seq <- Biostrings::readDNAStringSet(seq_path)

######################################################## MAIN CODE
### get data for plot 
# filter table
coverage_tb <- 
  coverage_tb_full %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "\\_bismark_bt2_pe.raw_coverage.bw")) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO")))

# create plot
baseplot <- 
  ggplot() +
  geom_rect(data = coverage_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = coverage), fill = "black") +
  facet_grid(rows = vars(genotype), 
             scales = "free_y") +
  ylab("") +
  xlab("") +
  scale_x_continuous(limits = c(0, width(cons_seq))) +
  # scale_y_continuous(limits = c(0, 90)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(plot = baseplot, 
       filename = file.path(outpath, str_c("LINE1.FLI_consensus.raw_coverage", "png", sep = ".")), 
       width = 15, height = 7.5)

