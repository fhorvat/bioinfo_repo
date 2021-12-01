a### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/methylation/Nov_2020/Bismark.unstranded_PE")

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
bw_path <- list.files(inpath, ".*\\.raw_coverage\\.bw$", full.names = T) 

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

######################################################## MAIN CODE
### get data for plot 
# filter table
coverage_tb <- 
  coverage_tb_full %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "\\.converted_bismark_bt2.bw|_bismark_bt2_pe.bw|_bismark_bt2_pe.raw_coverage.bw"), 
                genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"), 
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO"))) %>% 
  dplyr::arrange(genotype)

### plot for different repNames
purrr::map(unique(coverage_tb$repName), function(rmsk){
  
  # filter table 
  coverage_tb_repName <- 
    coverage_tb %>% 
    dplyr::filter(repName == rmsk) 
  
  # polygon total coverage
  pg <- 
    tibble(px = coverage_tb_repName$pos, 
           py = coverage_tb_repName$coverage, 
           sample_id = coverage_tb_repName$sample_id, 
           repName = coverage_tb_repName$repName, 
           genotype = coverage_tb_repName$genotype)
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_rect(data = pg, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "black") +
    facet_grid(rows = vars(genotype)) +
    ylab("") +
    xlab("") +
    # scale_x_continuous(limits = c(0, 6500)) +
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
         filename = file.path(outpath, str_c("raw_coverage", rmsk, "bw", "png", sep = ".")), 
         width = 15, height = 7.5)
  
  # return
  return(rmsk)
  
})
