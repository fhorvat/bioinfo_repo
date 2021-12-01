### INFO: 
### DATE: Mon Oct 29 16:48:51 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

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
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# miRBase gtf path
mirbase_gtf_path <- list.files(genome_path, pattern = "miRBase.*gff3", full.names = T)

# ensembl .gtf path
ensembl_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*gtf.gz", full.names = T)

# gene info path
gene_info_path <- list.files(path = genome_path, pattern = "ensembl.93.*[0-9]{6}.UCSC.*geneInfo.csv", full.names = T)

# miRBase mature miRNAs sequences path
mirbase_seq_path <- "/common/DB/genome_reference/miRBase/miRBase.22.20181029.mature.fa.gz"

######################################################## READ DATA
# read miRBase gff
mirna_gr <- rtracklayer::import.gff(con = mirbase_gtf_path)

# read ensembl .gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read gene info
gene_info <- readr::read_csv(gene_info_path)

# read miRBase sequences
mirbase_seq <- Biostrings::readRNAStringSet(filepath = mirbase_seq_path)

######################################################## MAIN CODE
##### 
## 3' UTR
# get 3' UTR coordinates, split by gene and reduce ranges
ensembl_3UTRs <-
  ensembl_gtf %>%
  gtfToGRanges(., filter = "three_prime_utr") %>% 
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F) %>%
  unlist(.)

# make names unique 
names(ensembl_3UTRs) <- make.unique(names(ensembl_3UTRs))

# get sequences of all 3' LTRs
ensembl_3UTRs_seq <- 
  BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, ensembl_3UTRs) %>% 
  as.character(.)

#####
## 7mers
kmer_7 <- 
  expand.grid(rep(list(c("A", "G", "T", "C")), 7)) %>% 
  do.call(str_c, .) %>% 
  DNAStringSet(.) %>% 
  reverseComplement(.) %>% 
  as.character(.)

# 
x <- purrr::map(names(ensembl_3UTRs_seq)[1:10], function(x) str_count(ensembl_3UTRs_seq[[x]], kmer_7))


lapply(names(ensembl_3UTRs_seq), function(x) str_count(ensembl_3UTRs_seq[[x]], kmer_rc))

kmer_7_df <- data.table(kmer = as.character(kmer_7), 
             kmer_rc = as.character(reverseComplement(kmer_7))) 

kmer_7_df[, names(ensembl_3UTRs_seq) := lapply(names(ensembl_3UTRs_seq), function(x) str_count(ensembl_3UTRs_seq[[x]], kmer_rc))]


# find 3' UTR matching nucleotides 2-7
match_7mer <- 
  Biostrings::vmatchPattern(pattern = DNAStringSet("TTTTTTT") %>% reverseComplement(.) %>% .[[1]], 
                            subject = ensembl_3UTRs_seq, 
                            max.mismatch = 0) %>% 
  unlist(.) %>% 
  length(.)

DNAString(kmer_7[1])

