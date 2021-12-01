### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS
# setdiff between two GRanges objects in a way which keeps names of final overlap
GRanges.parallelSetDiff <- function(gr1, gr2, ignore.strand = T){
  
  # get shared seqnames
  shared_seqnames <- 
    union(seqnames(gr1), seqnames(gr2)) %>% 
    as.character(.) %>% 
    unique(.)
  
  # change seqlevels both ranges share same seqlevels
  seqlevels(gr1) <- shared_seqnames
  seqlevels(gr2) <- shared_seqnames
  
  # find overlaps between two GRanges
  hits <- findOverlaps(gr1, gr2, ignore.strand = ignore.strand)
  
  # extract all overalaping features from subject as list
  grl <- extractList(gr2, as(hits, "List"))
  
  # parallel set difference - query vs. subject
  diff_list <- psetdiff(gr1, grl, ignore.strand = ignore.strand)
  
  # set names to diff. list, unlist
  names(diff_list) <- names(gr1)
  diff_gr <- unlist(diff_list)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# repeatMasker VIZ path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

# ensembl path
ensembl_path <- file.path(genome_dir, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.gtf.gz")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- read_delim(rmsk_path, delim = "\t")

# read gtf
ensembl_gtf <- read_delim(file = ensembl_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# clean repeatMasker, save
rmsk_filt <-
  rmsk_tb %>%
  dplyr::filter(repClass %in% c("LINE", "LTR")) %>%
  mutate(rmsk_id = str_c(seqnames, ":", start, "-", end)) %T>% 
  write_csv(., file.path(outpath, "rmsk.L1_and_LTRs.filtered.20190306.csv"))

# convert to GRanges
rmsk_gr <- 
  rmsk_filt %>%  
  GRanges(.)
names(rmsk_gr) <- rmsk_gr$rmsk_id

# get genes, reduce
genes_gr <- 
  ensembl_gtf %>% 
  gtfToGRanges(., filter = "gene") %>%
  GenomicRanges::reduce(., ignore.strand = T)

# remove whole genes, save
rmsk_without_genes <- GRanges.parallelSetDiff(gr1 = rmsk_gr, gr2 = genes_gr, ignore.strand = T) 

# save as .bed
rmsk_without_genes %>%
  sortSeqlevels(.) %>%
  sort(.) %T>%
  rtracklayer::export(object = ., con = file.path(outpath, "rmsk.L1_and_LTRs.filtered.20190306.removed_genes.bed"))

# add names as new column
mcols(rmsk_without_genes)$rmsk_id <- names(rmsk_without_genes)

# convert to table, join with details about each repeat
rmsk_filt_tb <- 
  rmsk_without_genes %>%
  as.tibble(.) %>% 
  left_join(., rmsk_filt %>% dplyr::select(rmsk_id, repName:repFamily), by = "rmsk_id")


### split to different elements, convert to GRanges
# LINE 1
line1_all <- 
  rmsk_filt_tb %>% 
  dplyr::filter(repFamily == "L1") %>% 
  dplyr::mutate(class = "LINE1_all") %>% 
  GRanges(.) 

# Line1 subfamilies = L1Md_T L1Md_A L1Md_Gf L1Md_F2 L1_Mus2 L1_Mus1 L1Md_F
line1_subfamilies <- 
  rmsk_filt_tb %>% 
  dplyr::filter(repName %in% c("L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F")) %>% 
  dplyr::mutate(class = repName) %>% 
  GRanges(.) 

# IAP all
iap_all <- 
  rmsk_filt_tb %>% 
  dplyr::filter(str_detect(repName, "IAP")) %>% 
  dplyr::mutate(class = "IAP_all") %>% 
  GRanges(.) 

# IAPEz-Int with LTRs IAPLTR1_Mm
iap_int_ltr <- 
  rmsk_filt_tb %>% 
  dplyr::filter(repName %in% c("IAPEy-int", "IAPLTR1_Mm")) %>% 
  dplyr::mutate(class = "IAPEy_int_IAPLTR1_Mm") %>% 
  GRanges(.)

# MuERV-int with MT2 LTRs = MERVL-int + MT2_Mm
mervl_mt2mm <- 
  rmsk_filt_tb %>% 
  dplyr::filter(repName %in% c("MERVL-int", "MT2_Mm")) %>% 
  dplyr::mutate(class = "MERVL_int_MT2_Mm") %>% 
  GRanges(.)

# MTA all LTRs and ints
mta_int_ltr <- 
  rmsk_filt_tb %>% 
  dplyr::filter(str_detect(repName, "MTA")) %>% 
  dplyr::mutate(class = "MTA_int_LTR") %>% 
  GRanges(.)

# ORR1A all ints and LTR subgroups
orr1a_int_ltr <- 
  rmsk_filt_tb %>% 
  dplyr::filter(str_detect(repName, "ORR1A")) %>% 
  dplyr::mutate(class = "ORR1A_int_LTR") %>% 
  GRanges(.)

### bind all and save
# bind rows to one table
retro_all_gr <- 
  c(line1_all, line1_subfamilies, 
    iap_all, iap_int_ltr, 
    mervl_mt2mm, mta_int_ltr, 
    orr1a_int_ltr) %T>%
  saveRDS(., file.path(outpath, "rmsk.L1_and_LTRs.filtered.20190306.removed_genes.GRanges.RDS"))

