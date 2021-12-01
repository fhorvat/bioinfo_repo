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

# classes of 21-23 nt alignments paths
class_21to23nt <- list.files(inpath, pattern = "smallRNA.21to23nt.*csv", full.names = T)

# classes of all lengths alignments paths
class_allLengths <- list.files(inpath, pattern = "smallRNA.all.*csv", full.names = T)

# histogram of alignment lengths
align_lengths <- list.files(file.path(inpath, "../alignment_length"), pattern = "align_length.*csv", full.names = T) 

# documentation path
documentation_path <- file.path(getwd(), "../../Data/Documentation")

# sample table path
sample_table_path <- list.files(path = documentation_path, ".*sampleTable.csv", full.names = T)

######################################################## READ DATA
# create empty table with all categories in order
category_df <- tibble(read_group = c("Mos_mRNA", "miRNA.sense", "protein_coding.sense", "rRNA", "SINE", "LINE", "LTR", "other_repeat", "annotated_pseudogene", "other", "not_annotated"))

# read library sizes, bind to one table
libsize_df <-
  purrr::map(align_lengths, function(lib_path) {
    
    # read library data.frame from files defined in libsize_list data.frame
    lib_df <-
      readr::read_csv(lib_path) %>%
      dplyr::mutate(ID = str_c(experiment, sample, sep = "."))
    
  }) %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(size_allLengths = sum(size),
                   size_21to23 = sum(size[read_length >= 21 & read_length <= 23])) %>%
  dplyr::mutate_at(vars(starts_with("size_")), funs(round(. / 1E6, 5)))

# read 21-23 nt counts, bind to one table
count_21to23 <-
  purrr::map(class_21to23nt, function(table_in){
    
    readr::read_csv(table_in) %>%
      dplyr::full_join(., category_df, by = "read_group") %>%
      dplyr::mutate(sample = replace(sample, is.na(sample), sample[1]),
                    experiment = replace(experiment, is.na(experiment), experiment[1]),
                    count = replace(count, is.na(count), 0)) %>% 
      dplyr::rename(count_21to23 = count)
    
  }) %>%
  dplyr::bind_rows(.) 

# read all size counts, bind to one table
count_allLengths <-
  purrr::map(class_allLengths, function(table_in){
    
    readr::read_csv(table_in) %>%
      dplyr::full_join(., category_df, "read_group") %>%
      dplyr::mutate(sample = replace(sample, is.na(sample), sample[1]),
                    experiment = replace(experiment, is.na(experiment), experiment[1]),
                    count = replace(count, is.na(count), 0)) %>% 
      dplyr::rename(count_allLengths = count)
    
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
  full_join(count_21to23, count_allLengths) %>%
  tidyr::unite(col = "ID", experiment, sample, sep = ".") %>%
  dplyr::left_join(., libsize_df, by = "ID") %>%
  dplyr::mutate(fpm_21to23 = count_21to23 / size_21to23, 
                fpm_allLengths = count_allLengths / size_allLengths, 
                ID = str_remove(ID, str_c(experiment_name, "."))) %>%
  dplyr::mutate_at(vars(starts_with("fpm")), funs(round(., 3))) %>% 
  dplyr::select(read_group, ID, count_21to23, count_allLengths, fpm_21to23, fpm_allLengths, size_21to23, size_allLengths)

# calculate average count/FPMs for whole experiment
count_avg <-
  count_df %>%
  dplyr::group_by(read_group) %>%
  dplyr::summarise(avg_count_allLengths = mean(count_allLengths) %>% round(., 2),
                   avg_count_21to23 = mean(count_21to23) %>% round(., 2),
                   avg_fpm_allLengths = mean(fpm_allLengths) %>% round(., 2),
                   avg_fpm_21to23 = mean(fpm_21to23) %>% round(., 2))

# gather, unite, spread to wide format
count_df_final <- 
  count_df %>%
  dplyr::select_at(vars(-starts_with("size"))) %>%
  tidyr::gather(variable, value, -(c(read_group, ID))) %>%
  tidyr::unite(temp, ID, variable) %>%
  tidyr::spread(temp, value) %>%
  dplyr::left_join(., count_avg, by = "read_group") %>%
  dplyr::right_join(., category_df, by = "read_group")

# write to separate sheets in xlsx 
write.xlsx(x = count_df_final %>% dplyr::select_at(vars("read_group", contains("count_all"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "count_all", 
           row.names = FALSE)
write.xlsx(count_df_final %>% dplyr::select_at(vars("read_group", contains("count_21to23"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "count_21to23", 
           append = TRUE, 
           row.names = FALSE)
write.xlsx(count_df_final %>% dplyr::select_at(vars("read_group", contains("fpm_all"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "fpm_all", 
           append = TRUE, 
           row.names = FALSE)
write.xlsx(count_df_final %>% dplyr::select_at(vars("read_group", contains("fpm_21to23"))) %>% as.data.frame(.), 
           file = file.path(outpath, str_c("smallRNA.", experiment_name, ".read_summary.xlsx")), 
           sheetName = "fpm_21to23", 
           append = TRUE, 
           row.names = FALSE)

# # calculate average count/FPMs for genotype/transfection
# count_avg_genotype <-
#   count_df %>%
#   dplyr::left_join(., sample_table, by = "ID") %>%
#   dplyr::group_by(read_group, strand, genotype, transfection) %>%
#   dplyr::summarise(avg_count_allLengths = mean(count_allLengths) %>% round(., 2),
#                    avg_count_21to23 = mean(count_21to23) %>% round(., 2),
#                    avg_fpm_allLengths = mean(fpm_allLengths) %>% round(., 2),
#                    avg_fpm_21to23 = mean(fpm_21to23) %>% round(., 2)) %T>%
#   readr::write_csv(x = ., path = file.path(outpath, str_c("smallRNA.", experiment_name, ".genotype_transfection.average.summary.csv")))



