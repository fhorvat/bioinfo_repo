### INFO: read smallRNA seq bam file, get counts over exons
### DATE: Wed May 16 02:54:16 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/class_reads")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(data.table)
library(xlsx)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "bamToGRangesList.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main outpath
outpath <- getwd()

# set main inpath
inpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# set experiment name
experiment_name <- "Eliska_mESC_MosIR"

# classes alignments paths
class_reads_path <- list.files(inpath, pattern = "smallRNA.*csv", full.names = T)

# library size
libsize_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/library_sizes/library_sizes.txt"

# documentation path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Documentation"

# sample table path
sample_table_path <- list.files(path = documentation_path, ".*sampleTable.csv", full.names = T)

######################################################## READ DATA
# create empty table with all categories in order
category_df <- tibble(read_group = c("Mos_mRNA", "miRNA.sense", "protein_coding.sense", "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", "other", "not_annotated"))

# read library sizes, bind to one table
libsize_df <- 
  readr::read_delim(file = libsize_path, delim = "\t", col_names = c("sample", "library_size")) %>% 
  dplyr::mutate(library_size = round(library_size / 1E6, 5))

# read class read counts, bind to one table
count_raw_df <-
  purrr::map(class_reads_path, function(table_in){
    
    readr::read_csv(table_in) %>%
      dplyr::full_join(., category_df, by = "read_group") %>%
      dplyr::mutate(sample = replace(sample, is.na(sample), sample[1]),
                    experiment = replace(experiment, is.na(experiment), experiment[1]),
                    count = replace(count, is.na(count), 0))
    
  }) %>%
  dplyr::bind_rows(.) 

# read sample table
sample_table <- 
  readr::read_csv(file = sample_table_path) %>% 
  dplyr::mutate(ID = str_remove(sample_id, ".SE")) %>% 
  dplyr::select(ID, genotype, transfection)

######################################################## MAIN CODE
# join tables with counts, join with library sizes, normalize
count_df <-
  count_raw_df %>% 
  dplyr::left_join(., libsize_df, by = "sample") %>%
  dplyr::mutate(fpm = (count / library_size) %>% round(., 3), 
                align_length = ifelse(str_detect(sample, "21to23nt"), "21to23nt", "all_lengths"), 
                sample = str_remove(sample, ".21to23nt")) %>%
  dplyr::select(read_group, sample, align_length, count, fpm, library_size)

# calculate average count/FPMs for whole experiment
count_avg <-
  count_df %>%
  dplyr::group_by(read_group, align_length) %>%
  dplyr::summarise(avg_count = mean(count) %>% round(., 3),
                   avg_fpm = mean(fpm) %>% round(., 3)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::gather(variable, value, -c(read_group, align_length)) %>% 
  tidyr::unite(temp, align_length, variable) %>% 
  tidyr::spread(temp, value)
  
# gather, unite, spread to wide format
count_df_final <- 
  count_df %>%
  dplyr::select(-library_size) %>%
  tidyr::gather(variable, value, -(c(read_group, sample, align_length))) %>%
  tidyr::unite(temp, sample, align_length, variable) %>%
  tidyr::spread(temp, value) %>%
  dplyr::left_join(., count_avg, by = "read_group") %>%
  dplyr::right_join(., category_df, by = "read_group")

# write to separate sheets in xlsx 
write.xlsx(x = count_df_final %>% dplyr::select_at(vars("read_group", contains("all_lengths_count"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "count_all", 
           row.names = FALSE)
write.xlsx(count_df_final %>% dplyr::select_at(vars("read_group", contains("21to23nt_count"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "count_21to23", 
           append = TRUE, 
           row.names = FALSE)
write.xlsx(count_df_final %>% dplyr::select_at(vars("read_group", contains("all_lengths_fpm"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "fpm_all", 
           append = TRUE, 
           row.names = FALSE)
write.xlsx(count_df_final %>% dplyr::select_at(vars("read_group", contains("21to23nt_fpm"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "fpm_21to23", 
           append = TRUE, 
           row.names = FALSE)



