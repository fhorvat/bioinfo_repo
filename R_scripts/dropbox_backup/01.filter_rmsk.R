### INFO: 
### DATE: Wed Apr 22 09:45:20 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome dir
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# repeatMasker path
rmsk_joined_path <- file.path(genome_dir, "rmsk.MesAur1.0.20160612.RefSeq.joined_rmsk_id.fa.out.gz")
rmsk_path <- file.path(genome_dir, "rmsk.MesAur1.0.20160612.RefSeq.clean.fa.out.gz")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read repeatMasker
rmsk_joined_tb <- readr::read_delim(rmsk_joined_path, delim = "\t")
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# get genes
genes_gr <- GRanges(genes_info)

# get all MT LTRs
MT_LTRs <- 
  rmsk_tb %>% 
  dplyr::filter(repClass == "LTR", str_detect(repName, "^MT"), !str_detect(repName, "-int")) %$% 
  repName %>% unique

# get MT2A
rmsk_joined_tb_filt <- 
  rmsk_joined_tb %>% 
  dplyr::filter(insertion_class %in% c("whole", "within")) %>% 
  dplyr::filter(str_detect(repName, "^MT2A.*MERVL_2A-int.*MT2A$")) %>% 
  dplyr::mutate(width = end - start + 1)

# # write
# rmsk_joined_tb_filt %>% 
#   dplyr::select(-insertion_class, -rmsk_id) %T>%
#   write_csv(., file.path(outpath, "rmsk.MesAur1.0.20160612.MT2A.MERVL_2Aint.csv"))


### find nearest gene
# rmsk GRanges
rmsk_gr <- GRanges(rmsk_joined_tb_filt)

# find nearest gene
genes_gr[GenomicRanges::nearest(rmsk_gr, genes_gr, ignore.strand = T)]

