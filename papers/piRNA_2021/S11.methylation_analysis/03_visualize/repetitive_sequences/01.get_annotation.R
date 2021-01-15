### INFO: 
### DATE: Wed Dec 09 13:55:37 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/methylation/repetitive")

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

library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# repeatMasker path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

# chosen LTR classes path
classes_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate"
classes_tb_path <- file.path(classes_tb_path, "all_LTR_classes 200730.xlsx")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path)

######################################################## MAIN CODE
# select relevant columns in LTR classes table
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

# get classes: Lx5, L1Mus, MuLV, MYSERVx, RLTR3x ERV1, IAP & ORR1
ltr_classes <- 
  classes_tb %>% 
  dplyr::filter(category_I %in% c("MuLV", "MYSERVx", "RLTR3x ERV1", "IAP", "ORR1"))
  
# get LINE1 classes
l1_classes <- 
  tibble(repName = c("Lx5", "L1_Mus3", "L1MdMus_I", "L1MdMus_II"), 
         repFamily = c("L1", "L1", "L1", "L1"), 
         category_I = c("Lx5", "L1Mus", "L1Mus", "L1Mus"), 
         type = c(NA, NA, NA, NA))
  
# get SINE classes
sine_classes <- 
  tibble(repName = c("B2_Mm1a", "B2_Mm1t", "B2_Mm2"), 
         repFamily = c("B2", "B2", "B2"), 
         category_I = c("SINEB2", "SINEB2", "SINEB2"), 
         type = c(NA, NA, NA))

# join to one table
all_classes <- rbind(ltr_classes, l1_classes, sine_classes)

# filter repeatMasker, get GRanges
rmsk_gr <- 
  rmsk_tb %>% 
  dplyr::select(-c(repClass, repFamily)) %>% 
  dplyr::right_join(., all_classes, by = "repName") %>% 
  GRanges(.)

# save as list of GRanges
split(rmsk_gr, rmsk_gr$category_I) %>% 
  saveRDS(., file = file.path(outpath, "rmsk.chosen_LTRs_and_L1s.GRanges.RDS"))
