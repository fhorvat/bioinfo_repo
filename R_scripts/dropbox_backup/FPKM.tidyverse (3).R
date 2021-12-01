### INFO: 
### DATE: Thu Oct 25 21:44:54 2018
### AUTHOR: Filip Horvat
# rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/other_projects/lab_meetings_Zagreb/2018_10_16")

######################################################## LIBRARIES
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)
library(tibble)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
se <- readRDS(file = "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.test.RDS")
colnames(se) <- str_remove(colnames(se), "\\.Aligned.sortedByCoord.out.bam")

sample_table <- readr::read_csv("sample_table.csv")
exons_table <- readr::read_csv("exons_width.csv")

######################################################## MAIN CODE
Rprof(memory.profiling = T)
### get FPKM from counts
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var = "gene_id") %>%
  tibble::as.tibble(.) %>%  
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>% 
  dplyr::left_join(sample_table %>% dplyr::select(-bam_path), by = "sample_id") %>%
  dplyr::left_join(exons_table, by = "gene_id") %>% 
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>%
  dplyr::select(gene_id, sample_id, fpkm, stage, genotype) %>% 
  tidyr::unite(stage_genotype, c("stage", "genotype"), sep = "_") %>% 
  dplyr::group_by(gene_id, stage_genotype) %>%
  dplyr::summarise(avg_fpkm = mean(fpkm) %>% round(., 3)) %>%
  dplyr::ungroup(.) %>% 
  tidyr::spread(key = stage_genotype, value = avg_fpkm)
Rprof(NULL)
x <- summaryRprof(memory = "both")$`by.self`
