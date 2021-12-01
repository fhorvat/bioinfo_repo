### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression")

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

library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# list of LINE1 full length elements path
line1_coords_path <- file.path(inpath, "Documentation", "L1s_nested_ours_20180517_PS_ZJM.xlsx")

######################################################## READ DATA
# read Zoe's list of LINE1 full length elements
line1_coords <- read.xlsx(line1_coords_path, sheet = "L1s_nested_ours_20180516_ZJM5")

######################################################## MAIN CODE
# tidy line1 data, save
line1_tidy <-
  line1_coords %>%
  magrittr::set_names(., make.unique(colnames(.))) %>%
  as.tibble(.) %>%
  dplyr::select(seqnames = chrName, start = chrStart, end = chrEnd, strand, id = uniqName, repClass, repFamily, repName) %>% 
  dplyr::mutate(strand = replace(strand, (is.na(strand) | strand == " "), "*")) %>% 
  distinct(.) %T>%
  write_csv(., file.path(outpath, "Documentation", "L1s_nested_ours_20180516.ZJM.tidy.csv"))

# save as bed
line1_gr <- 
  line1_tidy %>% 
  GRanges(.) %>% 
  sortSeqlevels(.) %>% 
  sort(.) %T>%
  rtracklayer::export(object = ., con = file.path(outpath, "Documentation", "L1s_nested_ours_20180516.ZJM.tidy.bed"))

# create genome GRanges
mm10_genome <- 
  tibble(seqnames = seqnames(BSgenome.Mmusculus.UCSC.mm10), 
         start = 1,
         end = seqlengths(BSgenome.Mmusculus.UCSC.mm10)) %>% 
  GRanges(.)

# get ranges not covered by LINE1s, save as bed
mm10_genome_masked <- 
  setdiff(mm10_genome, line1_gr, ignore.strand = T) %T>%
  rtracklayer::export(object = ., con = file.path(outpath, "Documentation", "L1s_nested_ours_20180516.ZJM.tidy.mm10_diff.bed"))

  
  