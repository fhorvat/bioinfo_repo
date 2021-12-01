### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Wed Nov 28 15:03:32 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis")

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

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR.new"

# library size path
library_size_path <- list.files(path = file.path(mapped_path, "4_library_size"), pattern = "library_hist.*mis_0.18to30nt.txt", full.names = T)

######################################################## READ DATA

######################################################## MAIN CODE
# read library size hist table, get different library sizes 
library_size_df <-
  purrr::map(library_size_path, function(library_path){
    
    readr::read_delim(library_path, delim = "\t", col_names = c("alignment_count", "alignment_length")) %>%
      dplyr::mutate(alignment_length = str_remove_all(alignment_length, "[0-9]*H|[0-9]*S|M")) %>%
      dplyr::group_by(alignment_length) %>%
      dplyr::summarise(alignment_count = sum(alignment_count)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(sample_id = basename(library_path) %>% str_remove_all(., "library_hist.|.txt"))
    
  }) %>%
  dplyr::bind_rows(.) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(library_size.mis_0.18to30nt = sum(alignment_count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE\\.mis_0\\.18to30nt")) %T>%
  readr::write_csv(., path = file.path(outpath, "library_size.Eliska_mESC_MosIR.csv"))

# ### write to other tables
# # IR expression
# write.xlsx(x = library_size_df %>%
#              as.data.frame(.),
#            file = file.path(outpath, "IR_expression", str_c("IR.", experiment, ".read_length.counts.xlsx")),
#            sheetName = "library_sizes",
#            append = TRUE,
#            row.names = FALSE)
# 
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




