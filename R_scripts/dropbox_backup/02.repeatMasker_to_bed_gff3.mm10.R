### INFO: 
### DATE: Sat Jul 20 21:51:54 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_locus/flanking_genes.mRNA/results/repeatMasker/mm10.inter_exons")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list all repeatMasker out files
rmsk_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/rmsk.mm10.20180919.raw.fa.out.gz"

# lnc1 in mm10 locus bed path
lnc1_bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_locus/flanking_genes.mRNA/results/mm10.inter_exons.bed"

######################################################## READ DATA
# read repeatMasker as table
rmsk_tb <- 
  readr::read_table2(file = rmsk_path, skip = 3, col_names = F) %>%
  dplyr::select(seqnames = X5, start = X6, end = X7, strand = X9, repName = X10, repClass_repFamily = X11, rmsk_id = X15) %>%
  tidyr::separate(col = repClass_repFamily, into = c("repClass", "repFamily"), sep = "/") %>%
  dplyr::mutate(strand = replace(strand, strand == "C", "-")) 

# read lnc1 in mm10 locus coordinates
lnc1_gr <- rtracklayer::import.bed(lnc1_bed_path)

######################################################## MAIN CODE
# read and tidy repeatMasker, convert to GRanges, subset to lnc1 locus
rmsk_gr <- 
  rmsk_tb %>% 
  GRanges(.) %>% 
  subsetByOverlaps(., lnc1_gr, ignore.strand = T)

# format and save as .bed
rmsk_gr %>% 
  as_tibble(.) %>% 
  dplyr::mutate(start = start - 1, name = repName, score = 1000, thickStart = start, thickEnd = end, itemRgb = "255,0,0") %>% 
  dplyr::select(chrom = seqnames, chromStart = start, chromEnd = end, name, score, strand, thickStart, thickEnd, itemRgb) %T>% 
  readr::write_delim(., file.path(outpath, "mm10.inter_exons.genomic.bed"), delim = "\t", col_names = F)

# save as gff3
rmsk_gr %>% 
  rtracklayer::export.gff3(., file.path(outpath, "mm10.inter_exons.genomic.gff3"))


