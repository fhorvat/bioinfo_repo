### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# list of experiments
experiment_list <- c("Deng_2014_Science_GSE45719", 
                     "Hamazaki_2015_Development_PRJDB2994", 
                     "Smallwood_2011_NatGenet_PRJEB2547", 
                     "Yamaguchi_2013_CellRes_GSE41908", 
                     "Gan_2013_NatCommun_GSE35005", 
                     "ENCODE_2014_Nature_GSE49417", 
                     "Veselovska_2015_GenomeBiol_GSE70116")

# accesory datasets path
accessory_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq"

# experiment paths
experiment_path <- file.path(accessory_path, experiment_list, "Data/Mapped/STAR_mm10")

# bam files
bam_path <- list.files(experiment_path, 
                       pattern = ".*\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam$|s_placenta_adult8wks_r.*PE\\.total\\.bam$|.*\\.bam$", 
                       full.names = T)

# add Fugaku bam files
bam_path_Fugaku <- list.files("/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Data/Mapped/STAR_mm10", 
                              pattern = ".*\\.bam$", 
                              full.names = T)

bam_path <- c(bam_path, bam_path_Fugaku)

######################################################## READ DATA

######################################################## MAIN CODE
# prepare sample table
sample_tb <- 
  tibble(bam_path = bam_path) %>% 
  mutate(sample_name = basename(bam_path),
         sample_id = str_remove(sample_name, "\\.genome\\.Aligned.*\\.bam$|\\.total\\.bam$|\\.bam"), 
         experiment = str_extract(bam_path, str_c(c(experiment_list, "Fugaku"), collapse = "|")), 
         stage = str_remove_all(sample_id, "^s_|_r[0-9]+\\.[S,P]E.*$")) %>% 
  dplyr::select(sample_id, stage, sample_name, experiment, bam_path) %T>%
  write_csv(., "developmental_stages.RNAseq.sample_table.20190501.csv")

