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
ensembl_3UTRs_seq <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, ensembl_3UTRs)

##### 
## miRNA
# prepare ranges of mature miRNA
mirna_gr <- mirna_gr[mcols(mirna_gr)$type == "miRNA"]
mcols(mirna_gr) <- mcols(mirna_gr)[, c("Name")]
names(mcols(mirna_gr)) <- "gene_id"

# get miRNA coordinates
mirna_coords <- mirna_gr[str_detect(mcols(mirna_gr)$gene_id, "mmu-let-7d-5p")][1]

# get let7a sequences
# mirna_seq <- BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, mirna_coords)
mirna_seq <- 
  mirbase_seq[str_detect(string = names(mirbase_seq), pattern = mcols(mirna_coords)$gene_id)] %>% 
  Biostrings::DNAStringSet(.) %>% 
  unlist(.)

# potential miRNA regulatory sites were found by searching the
# 3'UTR sequences for 7-nt matches, which included a 6-nt match to the miRNA seed
# (nucleotides 2-7) and either a seventh Watson-Crick match to miRNA nucleotide 8 or an
# adenosine opposite nucleotide 1

# find 3' UTR matching nucleotides 2-7
match_2to7 <- 
  Biostrings::vmatchPattern(pattern = subseq(mirna_seq, start = 2, end = 7) %>% reverseComplement(.), 
                            subject = ensembl_3UTRs_seq, 
                            max.mismatch = 0) %>% 
  unlist(.)


# # get gene info about matches
# match_info <- 
#   match_2to7 %>% 
#   as.tibble(.) %>% 
#   dplyr::mutate(names = str_remove(names, "\\..*$")) %>% 
#   # dplyr::distinct(names) %>% 
#   dplyr::left_join(., gene_info, by = c("names" = "gene_id")) %>% 
#   dplyr::filter(gene_name == "Hmga2")
# 
# ensembl_3UTRs_seq %>% .[names(.) == "ENSMUSG00000056758"] %>% subseq(., 59, 64)

# subset matching
ensembl_3UTRs_seq_subset <- ensembl_3UTRs_seq[names(ensembl_3UTRs_seq) %in% names(match_2to7)]
match_2to8 <- 
  Biostrings::vmatchPattern(pattern = subseq(mirna_seq, start = 2, end = 8) %>% reverseComplement(.), 
                            subject = ensembl_3UTRs_seq_subset, 
                            max.mismatch = 0) %>% 
  unlist(.)

subseq(ensembl_3UTRs_seq_subset, start = 1, end = 5)


