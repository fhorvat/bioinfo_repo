### INFO: 
### DATE: Fri Oct 09 09:41:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Raw/Cleaned")

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

library(networkD3)
library(htmlwidgets)
library(plotly)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Analysis")

# list trimming logs paths
trim_logs_path <- list.files(inpath, ".*\\.trim\\.log", full.names = T)

# mapping stats log
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi.trimmed"
mapping_log <- file.path(mapped_path, "bismark_summary_report.txt")

# deduplication stats logs
deduplication_log <- list.files(mapped_path, ".*\\.deduplication_report\\.txt", full.names = T)

######################################################## READ DATA
# read trimming logs
trim_logs_tb <- purrr::map(trim_logs_path, function(path){
  
  # read lines
  trim_log <- 
    readr::read_lines(path) %>%
    .[str_detect(., "in1=|Input:|KTrimmed:|Total Removed:")] %>% 
    tibble(raw = .) %>% 
    dplyr::mutate(clean = raw %>% str_remove_all(., ".*in1=|\t| reads.*|:| in2\\=.*|Total |,") %>% str_squish(.) %>% basename(.)) %>% 
    dplyr::select(clean) %>% 
    unique(.) %>% 
    dplyr::mutate(clean = ifelse(!str_detect(clean, " "), str_c("sample_id ", clean), clean)) %>% 
    tidyr::separate(clean, into = c("category", "value"), sep = " ") %>% 
    tidyr::pivot_wider(., names_from = "category", values_from = "value")
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "(?<=\\.[P,S]E).*")) %>% 
  dplyr::mutate_at(.vars = vars(!matches("sample_id")), .funs = as.numeric)

# read mapping logs
mapping_log_tb <- 
  readr::read_delim(mapping_log, delim = "\t") %>%
  dplyr::select(sample_id = File, pairs_into_mapping = `Total Reads`, mapped_reads = `Aligned Reads`) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "(?<=\\.[P,S]E).*"))

# read deduplication logs
deduplication_log_tb <- purrr::map(deduplication_log, function(path){
  
  # read and tidy
  dedup_log <- 
    readr::read_lines(path) %>% 
    .[str_detect(., "Total count of deduplicated leftover sequences|Total number of alignments analysed in")] %>% 
    tibble(raw = .) %>% 
    dplyr::mutate(value = str_remove(raw, ".*:") %>% 
                    str_squish(.) %>% 
                    str_remove_all(., " .*|\\t") %>% 
                    as.numeric(.),
                  category = c("dedup_in", "dedup_out"), 
                  sample_id = path %>% basename(.) %>% 
                    str_remove(., "_bismark_bt2_pe.deduplication_report.txt")) %>% 
    dplyr::select(sample_id, category, value) %>% 
    tidyr::pivot_wider(id_cols = sample_id, names_from = "category", values_from = "value")
  
}) %>%
  dplyr::bind_rows(.)

######################################################## MAIN CODE
# join tables 
bisulfite_stats <- 
  trim_logs_tb %>% 
  dplyr::left_join(., mapping_log_tb, by = "sample_id") %>%
  dplyr::left_join(., deduplication_log_tb, by = "sample_id") %>% 
  dplyr::mutate(trimmed_reads = Input -  Removed, 
                unmapped_reads = trimmed_reads - mapped_reads,
                removed_deduplication = dedup_in - dedup_out) %>% 
  dplyr::select(sample_id, 
                input_reads = Input, 
                removed_trimming = Removed, pass_trimming = trimmed_reads, 
                removed_mapping = unmapped_reads, pass_mapping = mapped_reads, 
                removed_deduplication, pass_deduplication = dedup_out) %T>% 
  readr::write_csv(., file.path(outpath, str_c("hamster_oocyte_Mov10l.bisulfite", 
                                               "no_conversion_filtering.library_stats.full.csv", 
                                               sep = ".")))

