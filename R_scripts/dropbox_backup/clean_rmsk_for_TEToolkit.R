### INFO: 
### DATE: Tue May 21 10:10:34 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Stein_2015_PLoSGenet_GSE57514/Analysis/TEToolkit_test")

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

library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# repeatMasker clean path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

# gtf downloaded from TE
gtf_TE_path <- file.path(inpath, "GRCm38_rmsk_TE.gtf.gz")
  
# UCSC seqnames path
UCSC_seqnames_path <- file.path(genome_dir, "GCF_000001635.20_GRCm38.ensembl2UCSC.txt")

######################################################## READ DATA
# # read repeatMasker
# rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read TE gtf
gtf_TE <- read_delim(file = gtf_TE_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))
  
# read UCSC seqnames
UCSC_seqnames <- read_delim(file = UCSC_seqnames_path, delim = "\t")

######################################################## MAIN CODE
# # save as gtf
# rmsk_tb_clean <- 
#   rmsk_tb %>% 
#   dplyr::filter(repClass %in% c("DNA", "DNA?", "LINE", "LINE?", "LTR", "LTR?", "Other", "RC", "RC?", "RNA", "Satellite", "SINE", "SINE?")) %>% 
#   dplyr::group_by(repName) %>%
#   dplyr::mutate(repName_count = 1:n()) %>% 
#   dplyr::ungroup(.) %>% 
#   dplyr::mutate(gene_id = repName,
#                 transcript_id = str_c(gene_id, "_dup", repName_count)) %>%
#   dplyr::select(seqnames:strand, gene_id, transcript_id, family_id = repFamily, class_id = repClass) %>%
#   dplyr::mutate(family_id = ifelse(is.na(family_id), class_id, family_id)) %>% 
#   GRanges(.)

# rtracklayer::export(object = rmsk_tb_clean, 
#                     con = file.path(outpath, "rmsk.mm10.20180919.gtf"),
#                     format = "gtf")

# change Ensembl seqnames to UCSC seqnames
ensembl_gtf_UCSC <-
  gtf_TE %>%
  dplyr::left_join(., UCSC_seqnames, by = c("X1" = "ensembl_name")) %>%
  dplyr::select(X1 = UCSC_name, X2:X9) %T>%
  write.table(x = ., file = file.path(outpath, "GRCm38_rmsk_TE.UCSC.gtf"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
