### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/summary.all")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# chosen LTR classes path
classes_tb_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/substitution_rate"
classes_tb_path <- file.path(classes_tb_path, "all_LTR_classes 200730.xlsx")

# expression path
expression_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression"
  
# RPM paths
oocyte_Mov10l1_path <- file.path(expression_path, "hamster_oocyte_Mov10l.RNAseq", "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.oocyte.RNAseq.RPM.csv")
# oocyte_CNOT6L_path <- file.path(expression_path, "hamster_oocyte_CNOT6L.RNAseq", "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.oocyte.RNAseq.RPM.csv")
testis_RNAseq_path <- file.path(expression_path, "hamster_testis_Mov10l.RNAseq", "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.testis.RNAseq.RPM.csv")
testis_small_RNAseq_path <- file.path(expression_path, "hamster_testis_Mov10l.small_RNAseq", "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.testis.small_RNAseq.sense.RPM.csv")

######################################################## READ DATA
# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(classes_tb_path)

# read RPMs
oocyte_Mov10l1 <- readr::read_csv(oocyte_Mov10l1_path)
# oocyte_CNOT6L <- readr::read_csv(oocyte_CNOT6L_path)
testis_RNAseq <- readr::read_csv(testis_RNAseq_path)
testis_small_RNAseq <- readr::read_csv(testis_small_RNAseq_path)

######################################################## MAIN CODE
### clean classes table
# select relevant columns
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

### clean expression data
# oocyte Mov10l1
oocyte_Mov10l1 %<>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Mov10l_WT", "Mov10l_KO"))) %>% 
  dplyr::arrange(genotype) %>% 
  tidyr::pivot_wider(., id_cols = repName, names_from = "genotype", values_from = "RPM", names_prefix = "oocyte.RNA.")

# # oocyte CNOT6L
# oocyte_CNOT6L %<>% 
#   tidyr::pivot_wider(., id_cols = repName, names_from = "genotype", values_from = "RPM", names_prefix = "oocyte.RNA.")

# testis RNA-seq
testis_RNAseq %<>% 
  dplyr::filter(age != "adult") %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_") %>% 
  dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_13dpp", "Mov10l_KO_13dpp", 
                                                               "Mov10l_WT_21dpp", "Mov10l_KO_21dpp"))) %>% 
  dplyr::arrange(genotype_age) %>% 
  tidyr::pivot_wider(., id_cols = repName, names_from = c("genotype_age"), values_from = "RPM", names_prefix = "testis.RNA.")

# testis small RNA-seq
testis_small_RNAseq %<>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_") %>% 
  dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_13dpp", "Mov10l_KO_13dpp", 
                                                               "Mov10l_WT_21dpp", "Mov10l_KO_21dpp"))) %>% 
  dplyr::filter(read_width >= 19, read_width <= 32) %>% 
  dplyr::group_by(genotype_age, repName) %>% 
  dplyr::summarise(RPM = sum(RPM)) %>% 
  dplyr::arrange(genotype_age) %>% 
  tidyr::pivot_wider(., id_cols = repName, names_from = c("genotype_age"), values_from = "RPM", names_prefix = "testis.sRNA.19to32nt.")

### add expression values to classes tb
# create one table
classes_tb_RPM <- 
  classes_tb %>% 
  dplyr::left_join(., oocyte_Mov10l1, by = "repName") %>% 
  dplyr::left_join(., testis_RNAseq, by = "repName") %>% 
  dplyr::left_join(., testis_small_RNAseq, by = "repName") %>% 
  dplyr::mutate_all(~(replace(., is.na(.), 0)))

# save
readr::write_csv(classes_tb_RPM, file.path(outpath, "all_LTR_classes.200730.all_RPMs.csv"))
  