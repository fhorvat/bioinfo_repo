### INFO: 
### DATE: Fri Jun 28 15:24:03 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("")

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

######################################################## SOURCE FILES

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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gtf path
gtf_path <- file.path(genome_dir, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.gtf.gz")

######################################################## READ DATA
# read .gtf
gtf_tb <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
# get coordinates of genes
genes_gr <- 
  gtfToGRanges(gtf_tb, filter = "gene") %>%
  GenomicRanges::split(., .$gene_id)

# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(gtf_tb, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = F, drop.empty.ranges = T)

# get coordinates of introns (= difference between full genes and exons)
introns_gr <- 
  GenomicRanges::setdiff(genes_gr, exons_gr, ignore.strand = F) %>% 
  GenomicRanges::reduce(., ignore.strand = F, drop.empty.ranges = T)

# remove empty ranges
introns_gr <- introns_gr[elementNROWS(introns_gr) > 0]


### remove all exons overlaping intronic regions, even from other genes
# get all exonic regions
exons_gr_all <- 
  exons_gr %>% 
  unlist(.) %>% 
  GenomicRanges::reduce(., ignore.strand = T, drop.empty.ranges = T)

# overlap introns with all exons
introns_gr_clean <- 
  introns_gr %>% 
  unlist(.) %>%  
  GRanges.parallelSetDiff(gr1 = ., gr2 = exons_gr_all, ignore.strand = T) %>% 
  GenomicRanges::split(., names(.))

# remove empty ranges
introns_gr_clean <- introns_gr_clean[elementNROWS(introns_gr_clean) > 0]

# save as .bed
introns_gr_clean %>%
  unname(.) %>% 
  unlist(.) %>% 
  sortSeqlevels(.) %>%
  sort(.) %T>%
  rtracklayer::export(object = ., con = file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.intronic.bed"))

 