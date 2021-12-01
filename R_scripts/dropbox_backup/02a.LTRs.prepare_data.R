### INFO: 
### DATE: Mon Oct 28 22:19:54 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/tmp/test")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# ensembl .gtf path
ensembl_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*[0-9]{6}.UCSCseqnames.gtf.gz"), full.names = T)

# rmsk path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# get only MT2_Mm
rmsk_mt2 <- 
  rmsk_tb %>% 
  dplyr::filter(repName %in% c("MT2A", "MT2B")) %>% 
  dplyr::select(chr = seqnames, ltr_start = start, ltr_end = end, ltr_strand = strand, 
                repName, repClass, repFamily, rmsk_id) %T>%
  readr::write_csv(., file.path(outpath, "rmsk.mm10.20180919.MT2A_B.fa.out.csv"))

# get and save exons on chromosome 1
reduced_exons <- 
  exons_gr %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(-width) %>% 
  dplyr::filter(seqnames == "chr1") %>% 
  dplyr::left_join(., genes_info %>% dplyr::select(gene_id, gene_biotype)) %T>%
  readr::write_csv(., file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chr1_reduced_exons.csv"))


### get unreduced exon coordinates
# get exon coordinates
exons_gr_all <- 
  ensembl_gtf  %>% 
  gtfToGRanges(., filter = "exon")

# get and save exons on chromosome 1
all_exons <- 
  exons_gr_all %>% 
  as_tibble(.) %>% 
  dplyr::select(seqnames:strand, gene_id, gene_biotype) %>% 
  dplyr::filter(seqnames == "chr1") %>% 
  readr::write_csv(., file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chr1_all_exons.csv"))

