### INFO: 
### DATE: Mon Aug 31 14:57:45 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/annotate_clusters")

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

library(openxlsx)
library(rtracklayer)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# read and tidy featureCounts counts
readFeatureCounts <- function(path){
  
  # read
  readr::read_delim(path, delim = "\t", comment = "#") %>% 
    set_colnames(., basename(colnames(.))) %>% 
    dplyr::select(-c(Chr:Length)) %>%
    dplyr::rename(gene_id = Geneid)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set experiment name
experiment_name <- "GL.mouse"

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gtf path
gtf_path <- file.path(genome_path, "ensembl.99.GRCm38.p6.20200415.UCSCseqnames.gtf.gz")

# gene info path
gene_info_path <- str_replace(gtf_path, "\\.gtf\\.gz$", ".geneInfo.csv")

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/GarciaLopez_2015_RNA_GSE59254"

# testis small RNA-seq expression path
testis_small_rnaseq_path <- file.path(base_path, "Analysis/expression.genes")
testis_small_rnaseq_sense_path <- list.files(testis_small_rnaseq_path, ".*\\.sense\\.counts\\.txt$", full.names = T)
testis_small_rnaseq_antisense_path <- list.files(testis_small_rnaseq_path, ".*\\.antisense\\.counts\\.txt$", full.names = T)

# sample table path
sample_table_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(sample_table_path, ".*\\.sampleTable\\.csv", full.names = T)

# library size path
library_size_path <- file.path(base_path, "Data/Mapped/STAR_mm10") 
library_size_path <- list.files(library_size_path, "library_sizes.txt", full.names = T, recursive = T)

######################################################## READ DATA
# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read small RNA-seq expression in testis
testis_smallrna_counts_sense <- readFeatureCounts(testis_small_rnaseq_sense_path)
testis_smallrna_counts_antisense <- readFeatureCounts(testis_small_rnaseq_antisense_path)

# read sample table
sample_tb <- readr::read_csv(sample_table_path)

# read library size table
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

######################################################## MAIN CODE
# get feature coordinates
features_tb <-
  testis_small_rnaseq_sense_path %>% 
  readr::read_delim(., delim = "\t", comment = "#") %>% 
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length)

# summarize library sizes
library_size_tb %<>% 
  dplyr::mutate(stage = str_extract(sample_id, "spermatozoa|spermatogonia|PGCs|MII|1C"), 
                read_length = str_extract(sample_id, "24to31nt|21to23nt|19to32nt"), 
                read_length = replace(read_length, is.na(read_length), "full")) %>% 
  dplyr::group_by(stage, read_length) %>% 
  dplyr::summarise(library_size = sum(library_size)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(read_length == "19to32nt") %>% 
  dplyr::mutate(sample_id = str_c("s_", stage)) %>% 
  dplyr::select(sample_id, library_size)

# transform counts to FPKM 
counts2FPM <- function(counts_tb){
  
  # long to wide
  counts_tb %>% 
    tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "count") %>% 
    dplyr::filter(str_detect(sample_id, "24to31nt")) %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE\\.24to31nt\\.bam$")) %>% 
    dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
    dplyr::left_join(., features_tb %>% dplyr::select(gene_id, width), by = "gene_id") %>% 
    dplyr::mutate(library_size = (library_size / 1E6), 
                  width = (width / 1E3), 
                  fpm = (count / library_size), 
                  fpkm = round((fpm / width), 4)) %>% 
    dplyr::select(gene_id, sample_id, fpm) %>% 
    tidyr::pivot_wider(., id_cols = gene_id, names_from = "sample_id", values_from = "fpm")
  
}

# get FPM
testis_smallrna_sense_fpm <- counts2FPM(testis_smallrna_counts_sense)
testis_smallrna_antisense_fpm <- counts2FPM(testis_smallrna_counts_antisense)

# save
readr::write_csv(testis_smallrna_sense_fpm, file.path(outpath, testis_small_rnaseq_sense_path %>% basename(.) %>% 
                                                        str_replace("\\.counts\\.txt", str_c(".", experiment_name, ".FPM.csv"))))
readr::write_csv(testis_smallrna_antisense_fpm, file.path(outpath, testis_small_rnaseq_antisense_path %>% basename(.) %>% 
                                                            str_replace("\\.counts\\.txt", str_c(".", experiment_name, ".FPM.csv"))))
