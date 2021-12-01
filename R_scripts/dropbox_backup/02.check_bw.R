#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: scales bigWig to RPMs
### DATE: Tue May 21 08:13:09 2019
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/datasets/Dicer_Mili_KO/Mapped/mm10_masked/1_mapped/test")

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
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(DESeq2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bw and bam path
bw_path <- file.path(inpath, "s_GV_DBL_old_r1.PE.perfect.bw")
bam_path <- file.path(inpath, "s_GV_DBL_old_r1.PE.perfect.bam")

######################################################## READ DATA
# read bigWig
bw <- rtracklayer::import(bw_path)

# read read counts
bam <- GenomicAlignments::readGAlignmentsList(file = bam_path, use.names = T, param = ScanBamParam(what = c("qname")))

######################################################## MAIN CODE
bam_coverage <- 
  bam %>% 
  coverage(.) %>% 
  GRanges(.)

# scale bigWig
mcols(bw)$score <- (mcols(bw)$score) / (round((read_stats$mapped_minus_rDNA / 1e6), 6))

# write bw
rtracklayer::export(object = bw, 
                    con = str_c(sample_id, ".scaled.bw"))


