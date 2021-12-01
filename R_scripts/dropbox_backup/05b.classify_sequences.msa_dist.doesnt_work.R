### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV/annotation")

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

library(BSgenome.Maur.UCSC.Siomi)
library(Biostrings)
library(msa)

# library(ape)
# library(kmer)
# library(dendextend)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# repeatModeler MYSERV annotation path
myserv_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV"
myserv_path <- file.path(myserv_path, "MYSERV.FLI_elements.bed")

######################################################## READ DATA
# read myserv coordinates
myserv_gr <- rtracklayer::import(myserv_path)

######################################################## MAIN CODE
# get sequences
insertions_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, myserv_gr)
names(insertions_seq) <- mcols(myserv_gr)$name

# do the MSA
seq_msa <- msa::msa(insertions_seq, method = "ClustalOmega")

# find distance
msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")
seq_dist <- seqinr::dist.alignment(msa_seqinr, "identity")

# build the neighbor-joining tree
insertions_phy <- ape::nj(seq_dist)

# ladderize tree
insertions_phy <- ape::ladderize(insertions_phy)

# convert to dendrogram
insertions_dnd <- stats::as.dendrogram(insertions_phy)

# plot
png(file = file.path(outpath, "MYSERV.FLI_elements.kmer6_distance.png"), width = 3000, height = 1000, units = "px")
plot(insertions_phy)
dev.off()

# cut the tree
insertions_clusters <- dendextend::cutree(insertions_phy, k = 3)
