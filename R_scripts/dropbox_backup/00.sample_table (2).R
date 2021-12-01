### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(tidyr)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# outpath
outpath <- getwd()

# inpath
inpath <- getwd()

# sample table
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation/lnc1_KO.RNAseq.20181211.sampleTable.raw.xlsx"

# raw files
raw_path <- list.files(path = "/common/RAW/Svoboda/lncRNA_KO/2018-12-11-CCRVDANXX", pattern = "*.txt.gz", full.names = T)

# barcodes path
barcodes_path <- "/common/RAW/Svoboda/lncRNA_KO/2018-12-11-CCRVDANXX/barcodes.txt"

# mapped files
bam_path <- list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Mapped/STAR_mm10", 
                       pattern = "*.genome.Aligned.sortedByCoord.out.bam$", full.names = T)
  
# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Mapped/STAR_mm10/log.Lnc1_KO.2018_Dec.stats_and_tracks.csv"
  
######################################################## READ DATA
# read sample table
sample_table <- openxlsx::read.xlsx(xlsxFile = sample_table_path, sheet = 2)

# # barcodes table
# barcode_table <- readr::read_delim(barcodes_path, delim = " ", col_names = c("sample_name", "freq", "barcode"))

# library size
library_size <- readr::read_csv(library_size_path)

######################################################## MAIN CODE
# set whether reads are single or paired end
pairing <- "SE"

# clean raw path
raw_df <- 
  tibble(raw_path) %>% 
  dplyr::mutate(sample_name = basename(raw_path) %>% str_remove(., "_sequence.txt.gz"), 
                raw_sample_id = str_remove(sample_name, "CCRVDANXX.*lane[4,5]"))

# clean mapped path
mapped_df <-  
  tibble(bam_path) %>% 
  mutate(sample_id = basename(bam_path) %>% str_remove_all(., ".genome.Aligned.sortedByCoord.out.bam"))
  
# clean sample table
sample_table_clean <-
  sample_table %>%
  as.tibble(.) %>%
  select(Sample.Name, barcode = Barcode) %>%
  mutate(genotype = str_extract(Sample.Name, "WT|Null"),
         stage = str_extract(Sample.Name, "GV|MII"),
         rep = str_extract(Sample.Name, "[0-9](?= )"),
         sample_id = str_c("s_Lnc1", stage, genotype, str_c("r", rep), sep = "_") %>% str_c(., ".SE")) %>%
  mutate(raw_sample_id = ifelse(test = str_detect(Sample.Name, "MII"),
                                yes = str_c("Lnc1MII", genotype, rep),
                                no = str_c("Lnc1GV", rep, genotype))) %>%
  left_join(., raw_df, by = "raw_sample_id") %>%
  left_join(., mapped_df, by = "sample_id") %>% 
  left_join(., library_size %>% select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id") %>% 
  select(sample_id, genotype, stage, rep, barcode, library_size, raw_sample_id, raw_path, bam_path) %T>%
  write_csv(., path = str_replace(sample_table_path, "raw\\.xlsx", "clean.csv")) %T>%
  write_csv(., path = file.path("/common/RAW/Svoboda/lncRNA_KO/2018-12-11-CCRVDANXX", str_replace(basename(sample_table_path), "raw\\.xlsx", "clean.csv")))


# # only for PE data
# if(pairing == "PE"){
#   
#   sample_table_clean %<>% 
#     rbind(., .) %>% 
#     dplyr::group_by(run) %>% 
#     dplyr::mutate(sequence = 1:n()) %>% 
#     dplyr::ungroup() %>% 
#     dplyr::mutate(run = str_c(run, "_", sequence), 
#                   sample_id = str_c(sample_id, "_", sequence)) %>% 
#     dplyr::select(-sequence)
#   
# }

# save make links script
sample_table_clean %>% 
  dplyr::arrange(sample_id) %>% 
  dplyr::mutate(sample_id = str_c(outpath, "/Data/Raw/Links/", sample_id, ".txt.gz")) %>% 
  dplyr::mutate(make_links = str_c("ln -s ", raw_path, " ", sample_id)) %$%
  make_links %T>% 
  readr::write_lines(., path = file.path(outpath, "Data/Raw/Links/make_links.sh"))
