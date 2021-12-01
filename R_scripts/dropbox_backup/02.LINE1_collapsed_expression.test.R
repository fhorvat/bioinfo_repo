### INFO: 
### DATE: Tue Jan 29 20:37:35 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# list of LINE1 full length elements path
line1_coords_path <- file.path(inpath, "Documentation", "LINE1_full_length.20180517.ZJM.tidy.csv")

# bam path
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Yang_2016_SciAdv_GSE83581/Data/Mapped/STAR_mm10.filtered/s_oocyte_r1.SE.19to32nt.bam" 

######################################################## READ DATA
# read Zoe's list of LINE1 full length elements
line1_coords <- read_csv(line1_coords_path)

######################################################## MAIN CODE
### prepare data
# form GRanges
line1_gr <- GRanges(line1_coords)

# tile each LINE1 to 50 tiles (= 2% of length)
line1_tiled <- tile(line1_gr, n = 50)
names(line1_tiled) <- line1_gr$id
line1_tiled <- unlist(line1_tiled)


### CNOT6L
# path
bam_path_CNOT6L <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/datasets/small_RNAseq/CNOT6L/Mapped/1_full_LINE1_reads/s_GV_WT_r1.PE.total.full_LINE1.bam"
bam_path_CNOT6L_unfiltered <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10/s_GV_WT_r1.PE.total.bam"

# read CNOT6L GV WT bam
bam_gr_CNOT6L <-
  readGAlignmentsList(file = bam_path_CNOT6L, use.names = TRUE)%>%
  unlist(.) %>% 
  sortSeqlevels(.) %>% 
  sort(.) %>% 
  grglist(.) %>%
  .[!duplicated(names(.))] %>%
  unlist(.)

# read CNOT6L GV WT bam
bam_gr_CNOT6L_unfiltered <-
  readGAlignmentsList(file = bam_path_CNOT6L_unfiltered, use.names = TRUE, 
                      param = ScanBamParam(which = reduce(line1_gr))) %>%
  unlist(.) %>% 
  sortSeqlevels(.) %>% 
  sort(.) %>% 
  grglist(.) %>%
  .[!duplicated(names(.))] %>%
  unlist(.)
