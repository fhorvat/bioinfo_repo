### INFO: 
### DATE: Mon Nov 05 12:55:11 2018
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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get chrY .gtf path
gtf_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chrYchrM.clean.gtf"

######################################################## READ DATA
# # read ensembl .gtf
# gtf_df <- 
#   read_delim(file = "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chrYchrM.gtf", 
#              delim = "\t", col_names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>% 
#   dplyr::mutate(attribute = str_remove(attribute, "exon_number \"\\d+\";")) %T>% 
#   write.table(x = ., file = file.path(dirname(gtf_path), "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.chrYchrM.clean.gtf"), 
#               col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# gtf_df <- 
#   read.table(file = gtf_path, sep = "\t", header = F) %>% 
#   as.tibble(.) %>% 
#   set_colnames(., c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

gtf_df <- read_delim(file = gtf_path, delim = "\t", 
                     col_names = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

######################################################## MAIN CODE
# Get the average width of last 2 exons on chromosome Y in each protein coding gene from xyz family. 
# If gene has only one exon, remove it from analysis. Be careful about strandness. 
gtf_clean <- 
  gtf_df %>% 
  dplyr::select(seqnames, feature, start, end, strand, attribute) %>% 
  dplyr::filter(feature == "exon", 
                seqnames == "chrY") %>% 
  dplyr::mutate(gene_id = str_extract(attribute, "(?<=gene_id \")ENSMUSG\\d+"), 
                gene_name = str_extract(attribute, "(?<=gene_name \")[:alnum:]+(?=\")"), 
                gene_biotype = str_extract(attribute, "(?<=gene_biotype \")[:alnum:]+_*[:alnum:]+(?=\")"))

gtf_clean$attribute %>% head