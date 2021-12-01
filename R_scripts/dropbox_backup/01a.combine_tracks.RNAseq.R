### INFO: combine track URLs and read stats table
### DATE: Mon Dec 14 00:09:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# get inpath
inpath <- getwd()

# get outpath
outpath <- getwd()

# list of experiments
experiment_list <- c("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq",
                     "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq",
                     "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.run_2.RNAseq")

# stats and tracks paths
tracks_path <- list.files(experiment_list, "log.tracks_URL.txt", full.names = T, recursive = T)
tracks_path <- tracks_path[str_detect(tracks_path, "STAR_Siomi.multimappers/5_merged_replicates")]

######################################################## READ DATA

######################################################## MAIN CODE
# read tracks
tracks_tbl <- purrr::map(tracks_path, function(path){
 
  # read table
  tracks_tb <- 
    readr::read_delim(file = path, col_names = "URL", delim = "\t") %>% 
    dplyr::mutate(perfect = str_detect(path, "perfect_reads"))
  
}) %>% 
  dplyr::bind_rows(.)

# clean table
tracks_tbl_tidy <-
  tracks_tbl %>% 
  dplyr::filter(str_detect(URL, "^http://hex.*")) %>% 
  dplyr::mutate(sample_id = basename(URL)) %>% 
  dplyr::filter(str_detect(sample_id, "^s_.*")) %>% 
  dplyr::mutate(experiment = URL %>% 
                  str_remove(., "http://hex.bioinfo.hr/~fhorvat/Svoboda/bw_tracks/in_house/hamster_KO/Siomi/") %>%
                  str_remove(., "/.*")) %>% 
  dplyr::mutate(perfect_2 = str_detect(URL, "perfect_reads"))

# check
if(any(tracks_tbl_tidy$perfect != tracks_tbl_tidy$perfect_2)) warning("Something's very wrong boss!")

# continue cleaning table
tracks_tbl_tidy_2 <- 
  tracks_tbl_tidy %>% 
  dplyr::select(-perfect_2) %>% 
  dplyr::filter(!str_detect(sample_id, "\\.bam\\.bai$")) %>% 
  dplyr::mutate(perfect = ifelse(perfect, "perfect", "mismatches"), 
                scaled = ifelse(str_detect(sample_id, "scaled"), "RPM_scaled", "raw"),
                file_type = ifelse(str_detect(sample_id, "\\.bw$"), "coverage", "individual_reads"), 
                sample_id = str_remove(sample_id, "\\.bam|\\.scaled\\.bw"),
                bw_name = str_c(str_remove(sample_id, "^s_"), perfect, sep = "."),
                bam_name = str_c(str_remove(sample_id, "^s_"), perfect, sep = "."), 
                URL = ifelse(test = (file_type == "coverage"),
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>% 
  dplyr::select(experiment, sample_id, scaled, file_type, perfect, URL)

# continue cleaning table
tracks_tbl_tidy_3 <- 
  tracks_tbl_tidy_2 %>% 
  tidyr::pivot_wider(id_cols = c(sample_id, experiment), names_from = c("scaled", "file_type", "perfect"), values_from = "URL") %>% 
  dplyr::select(experiment, sample_id, starts_with("RPM"), starts_with("raw"))

# save
readr::write_csv(tracks_tbl_tidy_3, file = file.path(outpath, "Mov10l1_hamster.all_multimappers.merged_replicates.RNAseq.tracks.20201214.csv"))








