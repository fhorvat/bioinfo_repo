### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/assemblies/golden_hamster.mesAur1/Mapped/CriGriPICR/LINE1/expand_10k.split/0_split_fasta/L1_1_CGr.1076274")

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

# deviate table path
bams_path <- list.files(inpath, ".*bam$", full.names = T) 

######################################################## READ DATA
# read data
coverage_tb_full <- purrr::map(bams_path, function(path){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = path, use.names = TRUE) %>% 
    unlist(.)
  
  # get coverage
  coverage_list <- 
    bam_gr %>% 
    coverage(.) %>% 
    as(., "IntegerList") %>% 
    unlist(.) %>% 
    split(., names(.)) 
  
  # set to table
  coverage_tb <- purrr::map(names(coverage_list), function(seqname){
    
    # create table
    cov_tb <- 
      coverage_list[[seqname]] %>% 
      tibble(coverage = .) %>% 
      dplyr::mutate(sample_id = basename(path), 
                    pos = 1:nrow(.)) %>% 
      dplyr::select(sample_id, pos, coverage)
    
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
  dplyr::mutate(sample_id = str_remove_all(sample_id, "\\.24to31nt\\.bam|\\.bam|s_mesAur_|.merged"))

# polygon total coverage
pg <- tibble(px = coverage_tb$pos, 
             py = coverage_tb$coverage, 
             sample_id = coverage_tb$sample_id) 

# create plot
baseplot <- 
  ggplot() +
  geom_rect(data = pg, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "black") +
  facet_grid(rows = vars(sample_id)) +
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
       filename = file.path(outpath, str_c("coverage", "png", sep = ".")), 
       width = 15, height = 7.5)

