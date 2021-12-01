### INFO: 
### DATE: Fri Dec 21 13:07:46 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Biostrings)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# differential expression analysis results path
results_gv_path <- file.path(inpath, "results", "Lnc1_KO.2018_Dec.diffExp.GRCm38.91.protein_coding.all.xlsx")
results_mii_path <- file.path(inpath, "results", "Lnc1_KO.2018_Dec.MII_Null_subset.diffExp.GRCm38.91.protein_coding.all.xlsx")

# ensembl .gtf path
ensembl_path <- file.path(genome_dir, "ensembl.91.GRCm38.p5.20180512.UCSCseqnames.gtf.gz")

######################################################## READ DATA
# read results
results_all <- 
  list(gv = read.xlsx(results_gv_path, sheet = "GV"), 
       mii_r1_r4 = read.xlsx(results_mii_path, sheet = "r1_r4"), 
       mii_r2_r3 = read.xlsx(results_mii_path, sheet = "r2_r3")) %>% 
  purrr::map(., as.tibble)

# read ensembl .gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# filter results, get all gene IDs
results_significant <- purrr::map(results_all, function(x) x %>% dplyr::filter(padj <= 0.1))
unique_gene_ids <- purrr::map(results_significant, function(x) x$gene_id) %>% unlist %>% unique

## 3' UTR
# get 3' UTR coordinates, split by gene and reduce ranges
ensembl_3UTRs <-
  ensembl_gtf %>%
  gtfToGRanges(., filter = "three_prime_utr") %>% 
  .[mcols(.)$gene_id %in% unique_gene_ids, ] %>% 
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F) %>%
  unlist(.)

# make names unique 
names(ensembl_3UTRs) <- make.unique(names(ensembl_3UTRs))

# get sequences of all 3' LTRs
ensembl_3UTRs_seq <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, ensembl_3UTRs)

### check for pattern
# get gene IDs of one stage
gene_ids_stage <- 
  results_significant[[1]] %$%
  gene_id

# get 3'UTRs
utrs_stage <- 
  ensembl_3UTRs_seq[str_remove(names(ensembl_3UTRs_seq), "\\.\\d+$") %in% results_stage] %>% 
  as.character(.)
  
# detect patterns in utrs
patterns <- 
  c("TTTTAAT", "TTTTATT", "TTTTATT") %>% 
  str_c(., collapse = "|")

# find patterns
str_count(string = utrs_stage, pattern = patterns)

# purrr::map(patterns, function(pattern) str_count(string = utrs_stage, pattern = pattern) %>% sum) 



