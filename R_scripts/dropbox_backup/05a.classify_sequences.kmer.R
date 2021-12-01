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
insertions_seq <- Biostrings::getSeq(BSgenome.Maur.UCSC.Siomi, myserv_gr)
names(insertions_seq) <- mcols(myserv_gr)$name

# convert to DNAbin
insertions_dnabin <- ape::as.DNAbin(insertions_seq)

# # cluster into OTUs
# insertions_otu <- otu(insertions_dnabin, k = 6, threshold = 0.97, method = "farthest", nstart = 20)

# calculate k-mer distance matrix
insertions_kdist <- kmer::kdistance(insertions_dnabin, k = 6)

# build the neighbor-joining tree
insertions_phy <- ape::nj(insertions_kdist)

# ladderize tree
insertions_phy <- ape::ladderize(insertions_phy)

# convert to dendrogram
insertions_dnd <- as.dendrogram(insertions_phy)

# plot
png(file = file.path(outpath, "MYSERV.FLI_elements.kmer6_distance.png"), width = 3000, height = 1000, units = "px")
plot(insertions_phy)
dev.off()

# cut the tree
insertions_clusters <- dendextend::cutree(insertions_dnd, k = 10)



### get the sequences of 1st cluster
# load additional libraries
library(msa)

### MSA each cluster
purrr::map(1:10, function(n){
  
  cat(n, "\n")
  
  # filter
  insertions_seq_filt <- insertions_seq[names(insertions_seq) %in% names(insertions_clusters[insertions_clusters == n])]
  
  # do the MSA
  seq_msa <- msa::msa(insertions_seq_filt, method = "ClustalOmega")
  
  # write as fasta
  seq_msa@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, str_c("MYSERV.ltr-1_family-16.cluster_", n, ".msa.fasta")))
  
  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)
  
  # save the consensus
  Biostrings::writeXStringSet(seq_consensus, file = file.path(outpath, str_c("MYSERV.ltr-1_family-16.cluster_", n, ".consensus.fasta")))
  
  
})
