### INFO: counts small RNA reads mapped to antisense features in repeatMasker 
### DATE: Tue Jul 31 23:18:59 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/repeat_expression.20180730/antisense_counts_rmsk")

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

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main outpath
outpath <- getwd()

# set main inpath
inpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# get repeatMasker path
rmsk_path <- list.files(path = genome_path, pattern = "rmsk.*[0-9]{6}.clean.fa.out.gz", full.names = T)

# list count vectors (.RDS)
count_list <- 
  list.files(path = inpath, pattern = "*.RDS", full.names = T) %>% 
  .[!str_detect(., "s_oocyte_19to24|s_oocyte_24to30")]

######################################################## READ DATA
# read repeatMasker
rmsk_df <- read_delim(file = rmsk_path, delim = "\t", col_types = cols(start = col_double(), end = col_double()))

######################################################## MAIN CODE
# tidy and filter repeatMasker (get only retrotransposons)
rmsk_df_filtered <- 
  rmsk_df %>% 
  tidyr::unite(col = coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::mutate(rmsk_idx = 1:nrow(.)) %>% 
  dplyr::filter(repClass %in% c("LINE", "LTR", "SINE"))

# read class read counts, bind to one table
count_df <-
  purrr::map(count_list, function(table_in){
    
    tibble(count = readRDS(table_in), 
           sample = basename(table_in) %>% str_remove_all(., "antisense.*(?=s_)|.RDS"), 
           # experiment = basename(table_in) %>% str_remove_all(., "antisense_count.rmsk.|.s_.*RDS$"),
           rmsk_idx = 1:nrow(rmsk_df)) %>% 
      dplyr::right_join(., rmsk_df_filtered %>% dplyr::select(rmsk_idx)) %>% 
      tidyr::separate(., col = sample, into = c("sample", "align_length"), sep = "\\.") %>% 
      dplyr::mutate(align_length = replace(align_length, is.na(align_length), "all_lengths")) %>% 
      dplyr::select(rmsk_idx, sample, align_length, count)
    
  }) %>%
  dplyr::bind_rows(.) 

# calculate average count/FPMs for whole experiment
count_avg <-
  count_df %>%
  dplyr::group_by(rmsk_idx, align_length) %>%
  dplyr::summarise(count = mean(count) %>% round(., 3)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(sample = "average") %>% 
  dplyr::select(rmsk_idx, sample, align_length, count)

# gather, unite, spread to wide format, join with repeatMasker
count_df_final <- 
  dplyr::bind_rows(count_df, count_avg) %>% 
  tidyr::unite(temp, sample, align_length) %>%
  tidyr::spread(temp, count) %>% 
  dplyr::filter(rowSums(.[, 2:ncol(.)]) > 0) %>% 
  left_join(., rmsk_df_filtered, by = "rmsk_idx") %>% 
  dplyr::select_at(vars(coordinates:repFamily, matches("s_MII.*|s_oocyte.*"), matches("average.*")))

# write separate read length tables
readr::write_csv(x = count_df_final %>%
                   dplyr::select_at(vars(coordinates:repFamily, matches(".*all_lengths"))) %>%
                   dplyr::arrange_at(vars(matches("average.*"))) %>%
                   dplyr::slice(nrow(.):(nrow(.) - 20000)), 
                 path = file.path(outpath, "antisense_reads.rmsk.all_lengths.counts_summary.csv"))

# write separate read length tables
readr::write_csv(x = count_df_final %>%
                   dplyr::select_at(vars(coordinates:repFamily, matches(".*21to23nt"))) %>%
                   dplyr::arrange_at(vars(matches("average.*"))) %>%
                   dplyr::slice(nrow(.):(nrow(.) - 20000)), 
                 path = file.path(outpath, "antisense_reads.rmsk.21to23nt.counts_summary.csv"))

# write separate read length tables
readr::write_csv(x = count_df_final %>%
                   dplyr::select_at(vars(coordinates:repFamily, matches(".*24to30nt"))) %>%
                   dplyr::arrange_at(vars(matches("average.*"))) %>%
                   dplyr::slice(nrow(.):(nrow(.) - 20000)), 
                 path = file.path(outpath, "antisense_reads.rmsk.24to30nt.counts_summary.csv"))


# # write to separate sheets in xlsx
# write.xlsx(x = count_df_final %>%
#              dplyr::select_at(vars(coordinates:repFamily, matches(".*all_lengths"))) %>%
#              dplyr::arrange_at(vars(matches("average.*"))) %>%
#              dplyr::slice(nrow(.):(nrow(.) - 20000)) %>%
#              as.data.frame(.),
#            file = file.path(outpath, "antisense_reads.rmsk.counts_summary.xlsx"),
#            sheetName = "all_length_reads_counts",
#            row.names = FALSE)
# 
# write.xlsx(x = count_df_final %>%
#              dplyr::select_at(vars(coordinates:repFamily, matches(".*21to23nt"))) %>%
#              dplyr::arrange_at(vars(matches("average.*"))) %>%
#              dplyr::slice(nrow(.):(nrow(.) - 20000)) %>%
#              as.data.frame(.),
#            file = file.path(outpath, "antisense_reads.rmsk.counts_summary.xlsx"),
#            sheetName = "21to23nt_reads_counts",
#            append = TRUE,
#            row.names = FALSE)
# 
# write.xlsx(x = count_df_final %>%
#              dplyr::select_at(vars(coordinates:repFamily, matches(".*24to30nt"))) %>%
#              dplyr::arrange_at(vars(matches("average.*"))) %>%
#              dplyr::slice(nrow(.):(nrow(.) - 20000)) %>%
#              as.data.frame(.),
#            file = file.path(outpath, "antisense_reads.rmsk.counts_summary.xlsx"),
#            sheetName = "24to30nt_reads_counts",
#            append = TRUE,
#            row.names = FALSE)



