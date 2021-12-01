### INFO: 
### DATE: Sat Apr 18 22:14:36 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/miRNA_targets_expression")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path(getwd(), "sequences")

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/pig/susScr11.Sscrofa11.1.GCA_000003025.6"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# get list of highest expressed transcript for each gene
transcripts_list_path <- file.path(inpath, "/transcripts_expression", "ensembl.93.Sscrofa11.1.20180920.UCSCseqnames.gene_transcript_expression.pig_GV.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read transcripts expression list
transcripts_list <- readr::read_csv(transcripts_list_path)

######################################################## MAIN CODE
# set animal
ensembl_name <- "sscrofa"

# Ensembl versions
ensembl_url <-
  tibble(ens_version = c(99, 98, 96, 95, 94, 93, 92, 91, 89, 86),
         date = c("Jan 2020", "Sep 2019", "Apr 2019", "Jan 2019", "Oct 2018", "Jul 2018", "Apr 2018", "Dec 2017", "May 2017", "Oct 2016"),
         URL_archive = c("www.ensembl.org",
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


### get sequences
# get 3'UTR sequences
UTR_seq <- 
  getSequence(id = transcripts_list$transcript_id, 
              type = "ensembl_transcript_id", 
              seqType = "3utr", 
              mart = mart) %>% 
  as_tibble(.) %>%
  dplyr::select(transcript_id = ensembl_transcript_id, seq_3UTR = `3utr`) %>%
  dplyr::filter(seq_3UTR != "Sequence unavailable") %>%
  dplyr::left_join(., transcripts_list, by = "transcript_id") %>%
  dplyr::select(gene_id, transcript_id, seq_3UTR)

# DNAStringSet
UTR_biostrings <- 
  UTR_seq$seq_3UTR %>% 
  Biostrings::DNAStringSet(.)
names(UTR_biostrings) <- str_c(UTR_seq$gene_id, UTR_seq$transcript_id, sep = ".")

# save fasta
Biostrings::writeXStringSet(x = UTR_biostrings, filepath = file.path(outpath, "ensembl.93.Sscrofa11.1.20200423.3UTR.fasta"))


