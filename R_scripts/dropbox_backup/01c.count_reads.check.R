### INFO: reproduces Rosa's Figure 4A from Genome Reasearch 2017 paper
### DATE: Mon Mar 11 15:47:01 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/Rosa")

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
library(GenomicFeatures)
library(GenomicAlignments)
library(Rsamtools)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# bam paths
# bam_path <- "%BAM_PATH"
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/Rosa/datasets/Gan_2013_NatCommun_GSE35005/filtered_reads/s_elongative_ST_r1.SE.bam"
bam_name <- basename(bam_path) %>% str_remove(., "\\.bam")

# range file path
range_path <- "/common/WORK/rosa/MT_elements/Data/MT_MT2_ORR_MLT_B1_B2_L1_IAP_elements.nogene.txt"

######################################################## READ DATA
# read in info for MT elements
mt <- read.delim(range_path, stringsAsFactors = FALSE)

# read bam
bam_gr <- readGAlignmentsList(bam_path, param = ScanBamParam(what = "qname"))

######################################################## MAIN CODE
# make ranges for MT elements
mtRanges <- GRanges(IRanges(mt$genoStart, mt$genoEnd), 
                    seqnames = Rle(mt$genoName), 
                    strand = mt$strand, 
                    name = mt$repName)

# split GRanges to GRangesList
rmsk_gr_list <- split(mtRanges, mtRanges$name)

# split GAlignments list to aligments of same read
bam_gr_split <- 
  bam_gr %>% 
  unlist(.) %>% 
  split(., mcols(.)$qname)

# count reads overlaping with repeatMasker = each read is counted only once for each of the classes in repeatMasker 
overlaps_table <- 
  findOverlaps(bam_gr_split, rmsk_gr_list, ignore.strand = T) %>% 
  as.data.table(.) %>% 
  .[, list(count = .N), by = subjectHits] %>% 
  as.tibble(.) %>% 
  dplyr::mutate(rmsk_id = names(rmsk_gr_list)[subjectHits], 
                sample_id = bam_name) %>% 
  dplyr::select(rmsk_id, count, sample_id) %>% 
  arrange(rmsk_id)

