### INFO: get expression in ENCODE mouse dataset
### DATE: Sun Jun 24 16:14:35 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)


######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS
# setdiff between two GRanges objects in a way which keeps names of final overlap
GRanges.parallelSetDiff <- function(gr1, gr2, ignore.strand = T){
  
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
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.91.*UCSCseqnames.reducedExons.RDS$", full.names = T)

# repeatMasker path
rmsk_path <- list.files(path = genome_dir, pattern = "rmsk.*.clean.fa.out.gz", full.names = T)

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read repeatMasker table
rmsk_df <-  readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# get GRanges from repeatMasker
rmsk_gr <- 
  rmsk_df %>% 
  GenomicRanges::makeGRangesFromDataFrame(.) %>% 
  GenomicRanges::reduce(.)

# substract ranges overlaping rmsk table
exons_filtered <- 
  exons_gr %>% 
  unlist(.) %>% 
  GRanges.parallelSetDiff(gr1 = ., gr2 = rmsk_gr, ignore.strand = T) %>% 
  GenomicRanges::split(., names(.))

# save filtered exons
saveRDS(object = exons_filtered, file = file.path(outpath, basename(exons_path) %>% str_replace(., ".RDS", ".rmskFiltered.RDS")))
