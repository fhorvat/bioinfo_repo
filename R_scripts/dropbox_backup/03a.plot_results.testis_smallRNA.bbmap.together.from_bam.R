### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/CriGriPICR.LINE1_consensus.testis_smallRNA.bbmap")

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

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/bbmap_mesAur1/4_library_size")

# sample table path
sample_tb_path <- list.files(documentation_path, pattern = "*.sampleTable.csv", full.names = T)

# stats and tracks path
stats_tb_path <- list.files(mapped_path, pattern = "library_sizes.txt", full.names = T)

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
      dplyr::mutate(repName = seqname, 
                    sample_id = basename(path), 
                    pos = 1:nrow(.)) %>% 
      dplyr::select(sample_id, repName, pos, coverage)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# read sample table
sample_tb <- readr::read_csv(sample_tb_path)

######################################################## MAIN CODE
# create table with genotype and library size
sample_tb_filt <-
  sample_tb %>% 
  dplyr::select(sample_id, genotype, age, library_size = count_19to32nt.Pepa) %>% 
  dplyr::mutate(library_size = (library_size / 10e5)) %>% 
  dplyr::mutate(genotype = str_remove(genotype, "Mov10l_"),
                genotype = factor(genotype, levels = c("WT", "KO"))) %>% 
  dplyr::filter(age == "13dpp")

# get unique repeat names
repeat_names <- unique(coverage_tb_full$repName)

### get data for plot 
purrr::map(repeat_names, function(repeat_name){
  
  # filter table
  coverage_tb <- 
    coverage_tb_full %>% 
    dplyr::filter(repName == repeat_name) %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "\\.24to31nt\\.bam")) %>% 
    dplyr::filter(sample_id %in% sample_tb_filt$sample_id)
  
  # limit L1-1_CGr family plot
  if(repeat_name == "L1-1_CGr"){
    coverage_tb %<>%
      dplyr::filter(pos < 6000)
  }
  
  # add genotype and library size to table, get average across genotype
  coverage_tb %<>% 
    dplyr::left_join(., sample_tb_filt, by = "sample_id") %>% 
    dplyr::select(sample_id, repName, pos, coverage, genotype, library_size) %>% 
    dplyr::mutate(coverage = round((coverage / library_size), 3)) %>%
    dplyr::group_by(genotype, pos) %>% 
    dplyr::summarise(coverage = mean(coverage)) %>% 
    dplyr::ungroup(.)
  
  # polygon total coverage
  pg <- tibble(px = coverage_tb$pos, 
               py = coverage_tb$coverage, 
               genotype = coverage_tb$genotype) 
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_rect(data = pg, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "lightgrey") +
    facet_grid(rows = vars(genotype)) +
    ylab("") +
    xlab("") +
    ggtitle(str_c(repeat_name, "24to31nt reads, testis 13dpp", sep = ", ")) +
    scale_x_continuous(limits = c(0, 6500)) +
    scale_y_continuous(limits = c(0, 10)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # save plot
  ggsave(plot = baseplot, 
         filename = file.path(outpath, str_c("LINE1", repeat_name, "CriGriPICR.5K_to_7k.ORFs.consensus", 
                                             "Mov10l1", "testis_smallRNA", "13dpp", "manual", "png", sep = ".")), 
         width = 15, height = 7.5)
  
  # return 
  return(pg)
  
}) %>% 
  bind_rows(.) %>% 
  dplyr::summarise(max_x = ceiling(max(px)),
                   max_y = ceiling(max(py)))
