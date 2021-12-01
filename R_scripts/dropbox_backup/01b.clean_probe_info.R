### INFO: 
### DATE: Fri Nov 02 19:59:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# set experiment path
exp_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Su_2004_ProcNatlAcadSciUSA_GSE1133"

# probe info path
probe_info_path <- file.path(exp_path, "GPL1073.GIN.gz")

######################################################## READ DATA
# read probe info, clean
probe_info <- 
  readr::read_delim(probe_info_path, delim = "\t", skip = 4) %>% 
  dplyr::select(index = Index, source = `Identifier Source`, probe_id = `Probe Set Name`, description = Description, 
                seq_type = `Sequence Type`, class1 = `Classification 1`, class2 = `Classification 2`) %>%
  dplyr::filter(seq_type == "mRNA") %>% 
  dplyr::select(-seq_type)

######################################################## MAIN CODE
# split probe data frame
probe_split <- split(probe_info, probe_info$source)

##### 
## names: Celera, GenBank, GNF, RefSeq, RTPS, UniGene

# Celera - not useful, no info about genes 
# probe_split[["GNF"]] %>% dplyr::filter_at(vars(starts_with("Classification")), any_vars(!is.na(.)))

# GenBank
genbank_df <- 
  probe_split[["GenBank"]] %>% 
  dplyr::mutate(gb = str_extract(description, "(?<=gb=)[:alnum:]+"), 
                gi = str_extract(description, "(?<=gi=)[:digit:]+"), 
                ug = str_extract(description, "(?<=ug=)Mm\\.[:digit:]+")) %>% 
  dplyr::select(-c(description, class1, class2))
  
# GNF
gnf_df <- 
  probe_split[["GNF"]]
  


