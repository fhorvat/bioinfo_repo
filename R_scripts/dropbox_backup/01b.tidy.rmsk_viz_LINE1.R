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

# repeatMasker VIZ path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

######################################################## READ DATA
# read repeatMasker
rmsk_tb <- read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# LINE1 families
line1_list <- c("L1_Mus1", "L1_Mus2", "L1Md_A", "L1Md_F", "L1Md_F2", "L1Md_T")

# tidy line1 data, save
line1_tidy <- 
  rmsk_tb %>% 
  dplyr::filter(repName %in% line1_list) %>% 
  mutate(seqnames = str_c(seqnames, ".", repName)) %>% 
  GRanges(.) %>% 
  sortSeqlevels(.) %>% 
  sort(.) %>% 
  reduce(., ignore.strand = F) %>% 
  as.tibble(.) %>% 
  separate(seqnames, c("seqnames", "repName"), "\\.") %>% 
  select(seqnames, start, end, strand, repName) %T>% 
  write_csv(., file.path(outpath, "Documentation", "LINE1_rmsk.20190213.tidy.csv"))

# remove repName from seqnames, save as .bed
line1_gr <- 
  line1_tidy %>% 
  GRanges(.) %T>% 
  rtracklayer::export(object = ., con = file.path(outpath, "Documentation", "LINE1_rmsk.20190213.tidy.bed"))

line1_gr %>% 
  as.tibble(.) %>% 
  group_by(repName) %>% 
  summarise(mean = mean(width), 
            sd = sd(width),
            min = min(width), 
            max = max(width))

# # create genome GRanges
# mm10_genome <-
#   tibble(seqnames = seqnames(BSgenome.Mmusculus.UCSC.mm10),
#          start = 1,
#          end = seqlengths(BSgenome.Mmusculus.UCSC.mm10)) %>%
#   GRanges(.)
# 
# # get ranges not covered by LINE1s, save as bed
# mm10_genome_masked <-
#   setdiff(mm10_genome, line1_gr, ignore.strand = T) %T>%
#   rtracklayer::export(object = ., con = file.path(outpath, "Documentation", "LINE1_full_length.20180517.ZJM.tidy.mm10_diff.bed"))


