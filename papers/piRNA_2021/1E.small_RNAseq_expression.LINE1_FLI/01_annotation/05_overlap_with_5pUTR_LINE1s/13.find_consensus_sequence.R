### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/05_overlap_with_5pUTR_LINE1s")

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

library(BSgenome.Maur.UCSC.Siomi)
library(seqinr)
library(Biostrings)
library(msa)
library(DECIPHER)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# insertions table path
insertions_path <- file.path(inpath, "LINE1.complete_ORFs.long_5pUTR.20200729.csv")

######################################################## READ DATA
# read insertion table
insertions_tb <- readr::read_csv(insertions_path)

######################################################## MAIN CODE
# to GRanges
insertions_gr <-
  insertions_tb %>%
  GRanges(.)
names(insertions_gr) <- str_c(mcols(insertions_gr)$repName, mcols(insertions_gr)$rmsk_id, sep = ".")
names(insertions_gr) <- str_replace_all(names(insertions_gr), "-", "_")

# get sequences
insertions_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, insertions_gr)

# # save as fasta
Biostrings::writeXStringSet(x = insertions_seq, filepath = file.path(outpath, str_c("LINE1.complete_ORFs.long_5pUTR.20200729", "fasta", sep = ".")))


### get consensus by repeat names
# split by repeat names
insertions_seq_list <- split(insertions_seq, insertions_gr$repName)

## MSA and get consesus
insertions_seq_consensus <- purrr::map(names(insertions_seq_list), function(seq_name){
  
  # do the MSA
  seq_msa <- msa::msa(insertions_seq_list[[seq_name]], method = "ClustalOmega")
  # seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq_list[[seq_name]], iterations = 100, refinements = 100, processors = 1, verbose = T)
  
  # write as fasta
  seq_msa@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, str_c("LINE1.complete_ORFs.long_5pUTR.20200729", seq_name, "msa", "fasta", sep = ".")))
  
  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)
  
  # return
  return(seq_consensus)
  
})

# join together
insertions_seq_consensus <- do.call(c, insertions_seq_consensus)
names(insertions_seq_consensus) <- names(insertions_seq_list)

# save the consensus
Biostrings::writeXStringSet(insertions_seq_consensus, file = file.path(outpath, str_c("LINE1.complete_ORFs.long_5pUTR.20200729", "consensus", "fasta", sep = ".")))


