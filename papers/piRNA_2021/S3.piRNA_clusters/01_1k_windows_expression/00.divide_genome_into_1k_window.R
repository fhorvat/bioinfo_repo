### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Analysis/expression.1k_window")

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

library(BSgenome.Maur.UCSC.MesAur1)
library(GenomicRanges)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

######################################################## READ DATA
# tile genome
genome_gr <- GenomicRanges::tileGenome(seqlengths = seqlengths(BSgenome.Maur.UCSC.MesAur1), 
                                       tilewidth = 1000, 
                                       cut.last.tile.in.chrom = T)
names(genome_gr) <- str_c(seqnames(genome_gr), ":", start(genome_gr), "-", end(genome_gr))

# export as .gtf
export.gff3(object = genome_gr, con = file.path(outpath, "MesAur1.1k_windows.gff"))

# # save as SAF
# line1_saf <-
#   line1_gr %>%
#   as_tibble(.) %>%
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>%
#   readr::write_delim(., file.path(outpath, "MesAur1.5k_windows.saf"), delim = "\t")

