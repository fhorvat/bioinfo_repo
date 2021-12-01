### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/Freimer_microarray/Sylamer_analysis")

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

library(GEOquery)
library(biomaRt)
library(beadarray)
library(beadarrayExampleData)
library(Biobase)
library(limma)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()


## results
# diff. expression results path
results_path <- file.path(inpath, "Freimer_2017.diff_exp.CrePos_vs_CreNeg.miR-15a.csv")


## gene annoation
# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# transcripts info path
transcripts_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.transcriptInfo.csv$"), full.names = T)

######################################################## READ DATA
# read results
results_tb <- readr::read_csv(results_path)

# read transcripts info
transcripts_info <- readr::read_csv(transcripts_info_path)

######################################################## MAIN CODE
### get all mouse 3' UTRs for Sylamer analysis
## load mart from Ensembl Biomart
# set animal
ensembl_name <- "mmusculus"

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


# get table of transcripts
utrs_tb_all <- 
  results_tb %>% 
  dplyr::filter(!is.na(gene_id)) %>% 
  left_join(., transcripts_info, by = "gene_id") 

# get 3'UTR sequences
utrs_seq_all <- getSequence(id = unique(utrs_tb_all$transcript_id), 
                            type = "ensembl_transcript_id", 
                            seqType = "3utr", 
                            mart = mart) 

# tidy 3'UTR sequences, get the longest 3'UTR per gene
utrs_tb_tidy <- 
  utrs_tb_all %>% 
  left_join(., utrs_seq_all, by = c("transcript_id" = "ensembl_transcript_id")) %>% 
  dplyr::rename(seq_3UTR = `3utr`) %>% 
  dplyr::filter(seq_3UTR != "Sequence unavailable") %>%
  dplyr::mutate(seq_length = nchar(seq_3UTR)) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::mutate(seq_length_rank = rank(-seq_length, ties.method = "random")) %>%
  dplyr::filter(seq_length_rank == 1) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(-seq_length_rank) %>% 
  dplyr::arrange(log2FoldChange)

# DNAStringSet
UTR_seq_biostrings <- 
  utrs_tb_tidy$seq_3UTR %>% 
  Biostrings::DNAStringSet(.)
names(UTR_seq_biostrings) <- utrs_tb_tidy$gene_name

# save fasta
Biostrings::writeXStringSet(x = UTR_seq_biostrings, 
                            filepath = file.path(outpath, "Freimer_2017.all_probes.ensembl_93.3pUTR.fasta"))

# save list of gene names ordered by logFC from upregulated to downregulated
utr_list <- 
  utrs_tb_tidy %$%
  gene_name %T>%
  readr::write_lines(., file.path(outpath, "Freimer_2017.all_probes.ensembl_93.3pUTR.logFC_ordered.txt"))

