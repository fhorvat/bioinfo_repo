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

# genome path
genome_path <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gtf path
gtf_path <- file.path(genome_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf")

# gene info path
gene_info_path <- str_replace(gtf_path, "\\.gtf$", ".geneInfo.csv")

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"

# testis small RNA-seq expression path
testis_small_rnaseq_path <- file.path(base_path, "Analysis/expression.genes")
testis_small_rnaseq_sense_path <- file.path(testis_small_rnaseq_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.sense.counts.txt")
testis_small_rnaseq_antisense_path <- file.path(testis_small_rnaseq_path, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.antisense.counts.txt")

# sample table path
sample_table_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- file.path(sample_table_path, "hamster_testis_Mov10l.20200528.sampleTable.csv")

# library size path
library_size_path <- file.path(base_path, "Data/Mapped/STAR_mesAur1.new/4_library_size") 
library_size_path <- file.path(library_size_path, "library_sizes.txt")
  
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
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l_KO|Mov10l_WT|Mov10l_HET"), 
                age = str_extract(sample_id, "13dpp|21dpp"), 
                read_length = str_extract(sample_id, "24to31nt|21to23nt|19to32nt"), 
                read_length = replace(read_length, is.na(read_length), "full")) %>% 
  dplyr::group_by(genotype, age, read_length) %>% 
  dplyr::summarise(library_size = sum(library_size)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(read_length == "19to32nt") %>% 
  dplyr::mutate(sample_id = str_c("s_testis_", genotype, "_", age)) %>% 
  dplyr::select(sample_id, library_size)

# transform counts to FPKM 
counts2FPM <- function(counts_tb){
  
  # long to wide
  counts_tb %>% 
    tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "count") %>% 
    dplyr::filter(str_detect(sample_id, "24to31nt")) %>% 
    dplyr::mutate(sample_id = str_remove(sample_id, "\\.24to31nt\\.bam$")) %>% 
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
readr::write_csv(testis_smallrna_sense_fpm, file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.sense.FPM.csv"))
readr::write_csv(testis_smallrna_antisense_fpm, file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.antisense.FPM.csv"))

