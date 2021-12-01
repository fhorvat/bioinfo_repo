### INFO: 
### DATE: Sun Dec 13 03:23:07 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Zhang_2021_NatCellBiol_GSE169528/Analysis/retrotransposon_expression")

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

# chosen LTR classes path
annotation_path <- file.path(inpath, "annotation")
annotation_path <- file.path(annotation_path, "all_LTR_classes 200730.xlsx")

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Zhang_2021_NatCellBiol_GSE169528"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi")
library_size_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(path = expression_path, pattern = ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read table with chosen LTR classes
classes_tb <- openxlsx::read.xlsx(annotation_path)

# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t") 

######################################################## MAIN CODE
### clean data
# select relevant columns in classes table
classes_tb %<>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repFamily, category_I, type) %>% 
  dplyr::mutate(repFamily = replace(repFamily, is.na(repFamily), "other"))

# filter sample table
sample_tb %<>% dplyr::select(sample_id, stage) 

# summarize library size table
library_size_tb %<>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id")


### get mean FPM expression
# read expression values
expression_tb_full <- purrr::map(expression_path, function(path){
  
  # get sample name
  sample_name <- path %>% basename(.) %>% str_remove(., "\\.expression\\.RDS")
  
  # read table
  tmp_tb <- 
    path %>% 
    readRDS(.) %>% 
    dplyr::mutate(sample_id = sample_name)
  
  # return
  return(tmp_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# add category to expression table, join with library size, get RPM
rpm_tb <- 
  expression_tb_full %>% 
  dplyr::left_join(., classes_tb, by = "repName") %>% 
  dplyr::group_by(sample_id, category_I) %>% 
  # dplyr::group_by(sample_id, repName) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"), 
                library_size = library_size / 1e6,
                RPM = count / library_size) %>%
  dplyr::select(sample_id, category_I, RPM)

# get mean RPM per age/genotype/experiment
expression_tb_mean <- 
  rpm_tb %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::mutate(stage = str_replace(stage, "PEA", "h post egg activation") %>% str_replace_all(., "_", " ")) %>% 
  dplyr::mutate(stage = factor(stage, levels = c("primary follicle", "secondary follicle", "GV", "MII",
                                                 "1C 9h post egg activation", 
                                                 "2C 33h post egg activation", "2C 44h post egg activation", "2C 52h post egg activation", 
                                                 "4C 44h post egg activation"))) %>% 
  dplyr::arrange(stage) %>% 
  dplyr::group_by(category_I, stage) %>% 
  dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = c(category_I), names_from = c("stage"), values_from = "RPM") %>% 
  dplyr::filter(!is.na(category_I)) %>% 
  dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0)))

### save 
# open workbook and sheets
wb <- createWorkbook("RPM")

# add sheet
addWorksheet(wb, sheetName = "RPM")

# write to workbook
writeDataTable(wb = wb, 
               sheet = "RPM", 
               x = expression_tb_mean)

# save workbook to disk
saveWorkbook(wb, file.path(outpath, str_c("all_LTR_classes.20210512", "long_RNAseq", "mean_RPM", "xlsx", sep = ".")), overwrite = TRUE)


