### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hubs/golden_hamster.Siomi/files/denovo_annotation")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# liftoff gtf path
gtf_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/Liftoff/MesAur1/RefSeq"
gtf_path <- file.path(gtf_path, "hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.gff")

######################################################## READ DATA
# read gtf using rtracklayer
gtf_gr <- rtracklayer::import.gff(gtf_path)

# # save as gff3
# rtracklayer::export.gff3(gtf_gr_clean, file.path(outpath, "hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.gff3"))

# # read gtf as table
# gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# remove begining from ID column
# gtf_gr_clean <- gtf_gr
# mcols(gtf_gr_clean)$ID <- str_remove(mcols(gtf_gr_clean)$ID, "^gene-|^rna-|^exon-|^cds-|^id-")
# mcols(gtf_gr_clean)$copy_id <- str_remove(mcols(gtf_gr_clean)$copy_id, "^gene-|^rna-|^exon-|^cds-|^id-")
# mcols(gtf_gr_clean)$Parent <- str_remove(mcols(gtf_gr_clean)$Parent, "^gene-|^rna-|^exon-|^cds-|^id-")
# 
# # remove tRNAs
# gtf_gr_clean <- gtf_gr_clean[mcols(gtf_gr_clean)$type != "tRNA"]


# get genes
gtf_genes <- 
  gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::select(seqnames, start, end, strand, ID, Dbxref, Name, gene, gene_biotype) %>% 
  dplyr::mutate(Dbxref = unlist(Dbxref))
