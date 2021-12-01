#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: scales bigWig to RPMs
### DATE: Tue May 21 08:13:09 2019
### AUTHOR: Filip Horvat

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Kabayama_2017_NucleicAcidsRes_PRJDB4628/Data/Mapped/STAR_mm10/4_merged_replicates/filter_by_length")

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

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# list of read stats, genome logs and merged logs
bw_path <- args$bw_path
read_stats_path <- args$read_stats_path

bam_path <- file.path(inpath, "s_Mili_WT.24to31nt.bam")
read_stats_path <- file.path(outpath, "..", "log.Kabayama_2017_NucleicAcidsRes_PRJDB4628.small_RNA.merged.stats_and_tracks.csv")

######################################################## READ DATA
# read bam
bam <- GenomicAlignments::readGAlignments(bam_path)

# read read counts
read_stats <- readr::read_csv(read_stats_path)

######################################################## MAIN CODE
# get sample name
sample_name <-
  basename(bam_path) %>%
  str_remove(., "(?<=SE).*\\.bam$|(?<=WT).*\\.bam$|(?<=KO).*\\.bam$")

# sum all the samples matching scaled one
scale_factor <-
  read_stats %>%
  dplyr::filter(str_detect(sample_id, sample_name)) %$%
  genome.mapped_minus_rDNA %>%
  sum(.)

# split by strand 
bam_strand <- list(plus = bam[strand(bam) == "+"], minus = bam[strand(bam) == "-"])

# get coverage, scale, save
bw_strand <- 
  purrr::map(c("plus", "minus"), function(strand){
  
  # get coverage
  bw <- 
    bam_strand[[strand]] %>% 
    coverage(.) %>% 
    GRanges(.)
  
  # add strand info
  strand(bw) <- ifelse(strand == "plus", "+", "-")
  
  # remove 0 score ranges
  bw <- bw[mcols(bw)$score != 0]
  
  # scale bigWig
  mcols(bw)$score <- (mcols(bw)$score) / (round((scale_factor / 1e6), 6))
  
  # negative score for minus strand
  if(strand == "minus"){
    mcols(bw)$score <- -mcols(bw)$score
  }

  # return 
  return(bw)
  
})

# join to one GRanges
bw_whole <- c(bw_strand[[1]], bw_strand[[2]])

# write bw
rtracklayer::export(object = bw_whole, con = str_c(sample_name, ".split_strand.scaled.bw"))
