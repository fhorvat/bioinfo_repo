### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/expression/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq")

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

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")

# expression files path
expression_path <- file.path(inpath, "expression_sum.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- 
  readr::read_delim(library_size_path, delim = "\t") %>% 
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA)

# read expression values
expression_tb <- purrr::map(expression_path, function(path){
  
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
  dplyr::bind_rows(.) %>% 
  dplyr::filter(sample_id != "s_testis_Mov10l1_WT_8.5dpp_So820-M12_r3.SE")

######################################################## MAIN CODE
### set element name
element_name <- "IAP"

### get expression 
# join with library size, get RPM
rpm_tb <- 
  expression_tb %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(library_size = library_size / 1e6,
                RPM = count / library_size) %>%
  dplyr::select(sample_id, RPM)

# join with sample table, summarize per genotype
rpm_mean <- 
  rpm_tb %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype) %>% 
  dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(genotype, "HET")) %>% 
  tidyr::pivot_wider(id_cols = c("genotype"), names_from = genotype, values_from = RPM, names_prefix = "RPM.") %>% 
  dplyr::mutate_all(.funs = ~(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(FLI_element = element_name) %>% 
  dplyr::select(FLI_element, everything())

# save
readr::write_csv(rpm_mean, file.path(outpath, str_c(element_name, "FLI_elements.RPM", "testis.8.5dpp", "run_2", "csv", sep = ".")))


