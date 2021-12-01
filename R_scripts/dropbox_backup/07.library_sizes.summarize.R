### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis")

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
### adds different library sizes (whole mapped library, 18 to 30 nt reads, 21 to 23 nt reads) to results tables
# experiment paths
experiment_paths <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
experiment <- experiment_list[1]

# get mapped path 
mapped_path <- experiment_paths[experiment]

# library size path
library_size_path <- list.files(path = file.path(mapped_path, "4_library_size"), pattern = "library_hist.*mis_[0-2]{1}.txt", full.names = T)

# read library size hist table, get different library sizes 
library_size_df <-
  purrr::map(library_size_path, function(library_path){
    
    readr::read_delim(library_path, delim = "\t", col_names = c("alignment_count", "alignment_length")) %>%
      dplyr::mutate(alignment_length = str_remove_all(alignment_length, "[0-9]*H|M")) %>%
      dplyr::group_by(alignment_length) %>%
      dplyr::summarise(alignment_count = sum(alignment_count)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(sample_id = basename(library_path) %>% str_remove_all(., "library_hist.|.txt"))
    
  }) %>%
  dplyr::bind_rows(.) %>% 
  tidyr::separate(sample_id, c("sample_id", "mismatch"), sep = ".mis_") %>% 
  dplyr::mutate(mismatch = str_c("mis_", mismatch)) %>% 
  tidyr::spread(mismatch, alignment_count) %>% 
  dplyr::mutate_at(vars(starts_with("mis")), funs(replace(., is.na(.), 0))) %>% 
  dplyr::mutate(mis_1 = mis_1 + mis_0, 
                mis_2 = mis_2 + mis_1) %>% 
  dplyr::arrange(sample_id, alignment_length) %>% 
  tidyr::gather(mismatch_allowed, alignment_count, mis_0:mis_2) %>% 
  dplyr::group_by(sample_id, mismatch_allowed) %>% 
  dplyr::summarise(sum_all = sum(alignment_count), 
                   sum_18to30nt = sum(alignment_count[(alignment_length >= 18) & (alignment_length <= 30)]), 
                   sum_21to23nt = sum(alignment_count[(alignment_length >= 21) & (alignment_length <= 23)])) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::gather(key = library_size, value = count, starts_with("sum_")) %>% 
  tidyr::unite(tmp, mismatch_allowed, library_size, sep = ".") %>% 
  tidyr::spread(tmp, count) %>% 
  dplyr::select(sample_id, str_c("mis_", rep(0:2, each = 3), ".sum_", c("all", "18to30nt", "21to23nt"))) %>% 
  data.table::setnames(x = ., old = colnames(.)[2:ncol(.)], new =  (str_remove(colnames(.)[2:ncol(.)], "sum_") %>% str_c(., "_reads"))) %T>%
  readr::write_csv(., path = file.path(outpath, "library_size", str_c("library_size.counts.", experiment, ".csv")))

### write to other tables
# IR expression
write.xlsx(x = library_size_df %>%
             as.data.frame(.),
           file = file.path(outpath, "IR_expression", str_c("IR.", experiment, ".read_length.counts.xlsx")),
           sheetName = "library_sizes",
           append = TRUE,
           row.names = FALSE)

# # hierarchically classified reads 
# write.xlsx(x = library_size_df %>%
#              as.data.frame(.),
#            file = file.path(outpath, "class_reads", str_c("class_reads.hierarchy.", experiment, ".whole_library.counts.xlsx")),
#            sheetName = "library_sizes",
#            append = TRUE,
#            row.names = FALSE)
# 
# # miRNA expression
# write.xlsx(x = library_size_df %>%
#              as.data.frame(.),
#            file = file.path(outpath, "miRNA_expression", str_c("mature_miRNA.", experiment, ".counts.xlsx")),
#            sheetName = "library_sizes",
#            append = TRUE,
#            row.names = FALSE)


                   
                            
