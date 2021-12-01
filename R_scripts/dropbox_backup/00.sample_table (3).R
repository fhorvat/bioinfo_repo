### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation")

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

# sample table path
sample_table_path_1 <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation/lnc1_KO.RNAseq.20181211.sampleTable.clean.csv"
sample_table_path_2 <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2017Sep_sequencing/Data/Documentation/RNAseq_2017_sampleTable.clean.csv"

# mapped files
bam_path_2 <- list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2017Sep_sequencing/Data/Mapped/STAR_mm10", 
                         pattern = "*.genome.Aligned.sortedByCoord.out.bam$", full.names = T)

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2017Sep_sequencing/Data/Mapped/STAR_mm10/log.2017Sep_sequencing.stats_and_tracks.csv"

######################################################## READ DATA
# read sample tables
sample_table_1 <- read_csv(sample_table_path_1)
sample_table_2 <- read_csv(sample_table_path_2)

# library size
library_size <- readr::read_csv(library_size_path)

######################################################## MAIN CODE
# clean mapped path
mapped_df <-  
  tibble(bam_path = bam_path_2) %>% 
  mutate(sample_id = basename(bam_path) %>% str_remove_all(., ".genome.Aligned.sortedByCoord.out.bam"))

# clean sample table 1
sample_table_1_clean <- 
  sample_table_1 %>% 
  mutate(experiment = "Lnc1_KO.2018_Dec") %>% 
  select(sample_id, experiment, genotype, stage, rep, barcode, library_size, raw_path, bam_path)

# clean sample table 2
sample_table_2_clean <- 
  sample_table_2 %>% 
  filter(sample_id %in% c("s_Lnc1Het_r1.SE", "s_Lnc1Het_r2.reseq.SE", "s_Lnc1Null_r1.SE", "s_Lnc1Null_r2.reseq.SE", "s_WT_r1.SE")) %>%
  left_join(., mapped_df, by = "sample_id") %>% 
  left_join(., library_size %>% select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id") %>% 
  mutate(genotype = str_remove(genotype, "Lnc1"), 
         rep = str_extract(sample_id, "[0-9](?=.reseq|.SE)"), 
         experiment = "Lnc1_KO.2017_Sep") %>% 
  select(sample_id, experiment, genotype, stage, rep, barcode, library_size, raw_path, bam_path)

# combine and save sample tables
rbind(sample_table_1_clean, sample_table_2_clean) %>% 
  write_csv(., file.path(outpath, "lnc1_KO.RNAseq.2017Sep_2018Dec.sampleTable.clean.csv"))
