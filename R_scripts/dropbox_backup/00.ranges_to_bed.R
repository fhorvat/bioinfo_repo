### INFO: 
### DATE: Tue Jan 22 18:30:44 2019
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
library(openxlsx)

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

# range file path
range_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/Rosa/MT_MT2_ORR_MLT_B1_B2_L1_IAP_elements.nogene.txt"

######################################################## READ DATA
# read in info for MT elements
mt <- read.delim(range_path, stringsAsFactors = FALSE)

######################################################## MAIN CODE
# make ranges for MT elements
mtRanges <- GRanges(IRanges(mt$genoStart, mt$genoEnd), 
                    seqnames = Rle(mt$genoName), 
                    strand = mt$strand, 
                    name = mt$repName)

# save as .bed
mtRanges %>%
  sortSeqlevels(.) %>%
  sort(.) %>% 
  reduce(., ignore.strand = T) %>% 
  .[start(.) > 0] %T>%
  rtracklayer::export(object = ., con = file.path(outpath, "MT_MT2_ORR_MLT_B1_B2_L1_IAP_elements.nogene.bed"))