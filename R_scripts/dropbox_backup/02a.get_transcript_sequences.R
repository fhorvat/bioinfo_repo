### INFO: gets 3'UTR sequences, longest for each gene
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/miR-205_pig")

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

library(biomaRt)
library(Biostrings)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 100

# genome path
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo.csv$"), full.names = T)

# transcripts info path
transcripts_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*\\.transcriptInfo\\.csv$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read transcripts info
transcripts_info <- readr::read_csv(transcripts_info_path)

######################################################## MAIN CODE
## load mart from Ensembl Biomart
# set animal
ensembl_name <- "sscrofa"

### download gene info from Ensembl
# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(100, 99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Apr 2020", "Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
                         "http://jan2020.archive.ensembl.org",
                         "http://sep2019.archive.ensembl.org",
                         "http://apr2019.archive.ensembl.org",
                         "http://jan2019.archive.ensembl.org",
                         "http://oct2018.archive.ensembl.org",
                         "http://jul2018.archive.ensembl.org",
                         "http://apr2018.archive.ensembl.org",
                         "http://dec2017.archive.ensembl.org",
                         "http://may2017.archive.ensembl.org",
                         "http://oct2016.archive.ensembl.org")) %>%
  dplyr::filter(ens_version == ensembl_version) %$%
  URL_archive

# load ENSEMBL mart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = stringr::str_c(ensembl_name, "_gene_ensembl"), host = ensembl_url)

# get 3'UTR sequences
utrs_seq_all <- 
  getSequence(id = genes_info$gene_id, 
              type = "ensembl_gene_id", 
              seqType = "3utr", 
              mart = mart) %>% 
  as_tibble(.) %>% 
  dplyr::rename(seq_3UTR = `3utr`) %>% 
  dplyr::filter(seq_3UTR != "Sequence unavailable") 

# get longest UTR for each gene
utr_seq_tidy <- 
  utrs_seq_all %>% 
  dplyr::mutate(width = nchar(seq_3UTR)) %>% 
  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::mutate(width_rank = rank(-width, ties.method = "random")) %>%
  dplyr::filter(width_rank == 1) %>% 
  dplyr::ungroup(.)

# DNAStringSet
UTR_seq_biostrings <- 
  utr_seq_tidy$seq_3UTR %>% 
  Biostrings::DNAStringSet(.)
names(UTR_seq_biostrings) <- utr_seq_tidy$ensembl_gene_id

# save fasta
Biostrings::writeXStringSet(x = UTR_seq_biostrings, 
                            filepath = file.path(outpath, str_c("ensembl", ensembl_version, "Sscrofa11.1.20200424.3pUTR.fasta", sep = ".")))

