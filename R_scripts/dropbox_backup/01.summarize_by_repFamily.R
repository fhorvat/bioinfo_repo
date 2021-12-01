### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/summary.repFamily")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# gene info path
gene_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression"
gene_info_path <- file.path(gene_info_path, "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.geneInfo.csv")

# expression table path
testis_tb_path <- file.path(inpath, "..", "hamster_testis_Mov10l.small_RNAseq", "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.testis.small_RNAseq.RPM.csv")
oocyte_tb_path <- file.path(inpath, "..", "hamster_oocyte_Mov10l.RNAseq", "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.oocyte.RNAseq.RPM.csv")

# choose testis or oocyte
expression_tb_path <- testis_tb_path

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read expression values
expression_tb <- readr::read_csv(expression_tb_path)

######################################################## MAIN CODE
### summarize gene info
# get repName and repClass combinations
gene_info_clean <- 
  gene_info %>% 
  dplyr::select(repName, repClass) %>% 
  unique(.)

### get expression 
# filter table by filter length, summarize
expression_sum <- 
  expression_tb %>% 
  dplyr::filter(read_width >= 24, read_width <= 31) %>% 
  dplyr::group_by(sample_id, repName) %>%
  dplyr::summarise(RPM = sum(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., gene_info_clean, by = "repName") %>% 
  dplyr::group_by(sample_id, repClass) %>%
  dplyr::summarise(RPM = sum(RPM)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(sample_id, "HET")) %>% 
  tidyr::pivot_wider(id_cols = c("sample_id"), names_from = c(repClass), values_from = RPM, names_prefix = "RPM.") %>% 
  dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))

# open workbook and sheets
wb <- createWorkbook("RPM")

# add sheet
addWorksheet(wb, sheetName = "RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "RPM", 
               x = rpm_filt)

# save workbook to disk
saveWorkbook(wb, file.path(outpath, str_c("rmsk_repFamily.200730", "small_RNA.testis", "24to31nt", "mean_RPM", "xlsx", sep = ".")), overwrite = TRUE)


