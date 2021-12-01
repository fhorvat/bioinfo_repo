### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/raw_rmsk.expression/hamster_testis_Mov10l.small_RNAseq")

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

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# expression files path
expression_path <- file.path(inpath, "expression_sum.sense.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

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
  dplyr::bind_rows(.)

######################################################## MAIN CODE
### get expression 
# add category to expression table, join with library size, get RPM
rpm_tb <- 
  expression_tb %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"), 
                library_size = library_size / 1e6,
                RPM = count / library_size) %>%
  dplyr::select(sample_id, repName, read_width, RPM) 

# join with sample table, summarize per genotype
rpm_mean <- 
  rpm_tb %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype, age, repName, read_width) %>% 
  dplyr::summarise(RPM = round(mean(RPM), 3)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(genotype, "HET")) 

# save
readr::write_csv(rpm_mean, "rmsk.Siomi.20200701.raw.SINE_LINE_LTRs.testis.small_RNAseq.sense.RPM.csv")
