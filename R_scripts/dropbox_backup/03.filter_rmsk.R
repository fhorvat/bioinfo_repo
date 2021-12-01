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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# repeatMasker VIZ path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# LINE1 row like in GenRes
# Line1 subfamilies = L1Md_T L1Md_A L1Md_Gf L1Md_F2 L1_Mus2 L1_Mus1 L1Md_F
# Lx subfamilies = Lx Lx2B2 Lx3_Mus Lx4B Lx5 Lx5c
# IAP all
# IAPEz-Int with LTRs IAPLTR1_Mm
# MuERV-int with MT2 LTRs = MERVL + MT2_Mm
# MTA all LTRs and ints
# ORR1A all ints and LTR subgroups

# tidy rmsk
rmsk_tidy <- 
  rmsk_tb %>% 
  dplyr::filter(repFamily %in% c("ERVL-MaLR", "L1", "ERVL", "ERVK")) %>%
  dplyr::filter(repFamily == "L1" | 
                  str_detect(repName, "IAP") | 
                  str_detect(repName, "MERVL|MT2_Mm") | 
                  str_detect(repName, "MTA") | 
                  str_detect(repName, "ORR1A")) %>% 
  dplyr::mutate(rmsk_id = str_c(seqnames, ":", start, "-", end)) %T>% 
  write_csv(., file.path(outpath, "rmsk.L1_and_LTRs.filtered.20190225.csv"))

