### INFO: reads gtf, gets all exons of one gene and reduces ranges, saves object as RDS
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/chinese_hamster/CHOK1GS_HDv1.GCA_900186095.1")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "gtfToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# gtf path
gtf_path <- list.files(path = inpath, pattern = "ensembl.*UCSCseqnames.gtf.gz")
# gtf_path <- list.files(path = inpath, pattern = "refGene.*.gtf.gz")
gtf_path <- list.files(path = inpath, pattern = "ensembl.91.*.gtf.gz")

######################################################## READ DATA
# read gtf
gtf_df <- read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
### for each gene get reduced exonic regions
# convert GTF to GRanges, get only exons, reduce
exons_gr <-
  gtfToGRanges(gtf_df, filter = "exon") %>%
  GenomicRanges::split(., .$gene_id) %>%
  GenomicRanges::reduce(., ignore.strand = T)

# save RDS
saveRDS(object = exons_gr, file = file.path(outpath, stringr::str_replace(basename(gtf_path), ".gtf.gz", ".reducedExons.RDS")))
