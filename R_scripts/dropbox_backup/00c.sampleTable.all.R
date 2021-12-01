### INFO: get expression of all genes in lncKO data
### DATE: Wed May 23 19:20:31 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/diffExp")
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set input path
inpath <- getwd()

# set output path
outpath <- getwd()

# get gene info path
info_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.geneInfo.csv", full.names = T)

# get paths of reduced exons
exons_path <- list.files(genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS", full.names = T)

# stats and tracks, bam paths
mapped_path <- c("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2016Nov_sequencing/Data/Mapped/STAR_mm10", 
                 "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2017Sep_sequencing/Data/Mapped/STAR_mm10")

# sample table path
sample_table_path <- c("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2016Nov_sequencing/Data/Documentation/RNAseq_2016_11_23_sampleTable.clean.csv", 
                       "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/2017Sep_sequencing/Data/Documentation/RNAseq_2017_sampleTable.clean.csv")

######################################################## READ DATA
# read stats and tracks
stats_df <- 
  list.files(mapped_path, pattern = "log.*stats_and_tracks.csv", full.names = T) %>% 
  purrr::map(., readr::read_csv) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(experiment = str_remove(experiment, "_sequencing"), 
                sample_id = str_c(sample_id, ".", experiment))

# read sample tables
sample_df <- 
  purrr::map(sample_table_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = str_c(sample_id, ".", experiment))

# create bam files table
bam_df <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>% 
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".genome.Aligned.sortedByCoord.out.bam", ""), 
                experiment = stringr::str_extract(bam_path, "2016Nov|2017Sep")) %>% 
  dplyr::mutate(sample_id = str_c(sample_id, ".", experiment))

######################################################## MAIN CODE
# construct sample table
sample_table <-
  sample_df %>% 
  dplyr::select(sample_id, genotype, resequencing) %>% 
  dplyr::left_join(., bam_df, by = "sample_id") %>% 
  dplyr::left_join(., stats_df %>% dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA), by = "sample_id") %>% 
  dplyr::mutate(genotype = replace(genotype, genotype == "Lnc12", "Lnc12Null"), 
                genotype = replace(genotype, genotype == "Lnc21", "Lnc21Null")) %>% 
  dplyr::select(sample_id, genotype, library_size, bam_path) %T>%
  readr::write_csv(., path = file.path(outpath, "lncKO.all.sample_table.csv"))

