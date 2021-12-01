### INFO: read runInfo.txt from SRA and creates renaming script for fastq files
### DATE: 28. 11. 2017.  
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/Papd7_KO/retrotransposon_expression_signature")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# inpath
inpath <- getwd()

# outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# repeat masker path
rmsk_path <- list.files(genome_path, "rmsk.*\\.clean\\.fa\\.out\\.gz", full.names = T)

# joined rmsk ID path
rmsk_joined_path <- list.files(genome_path, "rmsk.*\\.joined_rmsk_id\\.fa\\.out\\.gz", full.names = T)

# gtf path
gtf_path <- list.files(genome_path, "ensembl\\.99.*\\.UCSCseqnames.geneInfo.csv", full.names = T)

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read joined repeatMasker
rmsk_joined <- readr::read_delim(rmsk_joined_path, delim = "\t")

# read gene info
gene_info <- readr::read_csv(gtf_path)

######################################################## MAIN CODE
# filter repeatMasker
rmsk_tb_filt <- 
  rmsk_tb %>% 
  dplyr::filter(repClass %in% c("LINE", "SINE", "LTR", "Satellite")) %>% 
  dplyr::mutate(gene_id = str_c("rmsk_", rmsk_ID)) %>%
  dplyr::mutate(type = "exon", 
                insertion_class = "single") %>% 
  dplyr::rename(rmsk_id = rmsk_ID)
  
# filter joined repeatMasker
rmsk_joined_filt <- 
  rmsk_joined %>% 
  dplyr::filter(rmsk_id %in% rmsk_tb_filt$rmsk_id) %>% 
  dplyr::mutate(gene_id = str_c("rmsk_", rmsk_id)) %>%
  dplyr::mutate(type = "gene")

# join, create GRanges
rmsk_all <- 
  dplyr::bind_rows(rmsk_tb_filt, rmsk_joined_filt) %>% 
  dplyr::mutate(strand = replace(strand, strand %in% c("-/+", "+/-"), "*")) %>% 
  dplyr::arrange(rmsk_id) %>% 
  GRanges(.)

# write as .gtf
rtracklayer::export.gff(rmsk_all, file.path(genome_path, "rmsk.mm10.20180919.joined_rmsk_id.gtf"))

# get gene info
rmsk_gene_info <- 
  rmsk_all %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::mutate(gene_name = gene_id, 
                gene_biotype = repName, 
                gene_description = str_c(gene_id, insertion_class, sep = ",")) %>% 
  dplyr::select(gene_id, seqnames, start, end, strand, gene_name, gene_biotype, gene_description)

# save
readr::write_csv(rmsk_gene_info, file.path(genome_path, "rmsk.mm10.20180919.joined_rmsk_id.geneInfo.csv"))

