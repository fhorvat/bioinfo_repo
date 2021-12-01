### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/class_reads")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(xlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA

######################################################## MAIN CODE
### summarize read classes for all samples in experiment 
# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.all_plasmids") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  experiment <- "NIH3T3_transfected.2018"
  
  # get mapped path 
  mapped_path <- experiment_paths[experiment]
  
  # classified reads path
  read_class_path <- list.files(path = file.path(mapped_path, "6_filtered_21to23nt"), pattern = "rmsk.read_class.*RDS", full.names = T)
  
  # read class reads table, gets different mismatch counts
  read_class_df <-
    purrr::map(read_class_path, function(x){
      
      readRDS(x)[["repFamily"]] %>% 
        dplyr::filter(!(repFamily == "not_repeat"))
      
    }) %>%
    dplyr::bind_rows(.) %>% 
    tidyr::spread(key = sample_id, value = count) %>% 
    dplyr::mutate_at(vars(matches("^s_.*")), funs(replace(., is.na(.), 0)))
  
  # write to separate sheets in xlsx 
  write.xlsx(x = read_class_df %>%
               as.data.frame(.), 
             file = file.path(outpath, str_c("rmsk.class_reads.", experiment, ".21to23nt.counts.xlsx")), 
             sheetName = "allowed_0_mismatches", 
             row.names = FALSE)
  
}

