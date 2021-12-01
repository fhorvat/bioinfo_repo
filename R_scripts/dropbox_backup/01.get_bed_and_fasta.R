### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/mouse.mm10/FLI_elements/LINE1/RNAi_paper_list")

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

library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)
library(msa)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# LINE1s rmsk path
line1_path <- file.path(inpath, "LINE1.full_length.RNAi_2019_paper.csv")

######################################################## READ DATA
# read LINE1 table
line1_tb <- readr::read_csv(line1_path)

######################################################## MAIN CODE
# to GRanges
line1_gr <- 
  line1_tb %>% 
  dplyr::select(seqnames = chrName, start = chrStart, end = chrEnd, strand = strand, repName) %>% 
  dplyr::mutate(rmsk_id = str_c(repName, "|", seqnames, ":", start, "-", end)) %>% 
  dplyr::filter(!is.na(strand)) %>% 
  GRanges(.)
names(line1_gr) <- mcols(line1_gr)$rmsk_id

# save as .bed
rtracklayer::export.bed(line1_gr, file.path(outpath, "LINE1.full_length.RNAi_2019_paper.bed"))

# # save as SAF
# line1_saf <- 
#   line1_gr %>% 
#   as_tibble(.) %>% 
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
#   readr::write_delim(., file.path(outpath, "LINE1.full_length.RNAi_2019_paper.saf"), delim = "\t")

# get sequences
line1_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, line1_gr)

# save as fasta
Biostrings::writeXStringSet(x = line1_seq, filepath = file.path(outpath, "LINE1.full_length.RNAi_2019_paper.fasta"))


### get consensus by repeat names
# split by repeat names
line1_seq_list <- split(line1_seq, line1_gr$repName)

## MSA and get consesus
line1_seq_consensus <- purrr::map(names(line1_seq_list), function(seq_name){
  
  # do the MSA
  seq_msa <- msa::msa(line1_seq_list[[seq_name]], method = "ClustalOmega")
  # seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq_list[[seq_name]], iterations = 100, refinements = 100, processors = 1, verbose = T)
  
  # # write as fasta
  # seq_msa@unmasked %>%
  #   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("LINE1.full_length.RNAi_2019_paper", seq_name, "msa.fasta", sep = ".")))
  
  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)
  
  # return
  return(seq_consensus)
  
})

# join together
line1_seq_consensus <- do.call(c, line1_seq_consensus)
names(line1_seq_consensus) <- names(line1_seq_list)

# save the consensus
Biostrings::writeXStringSet(line1_seq_consensus, file = file.path(outpath, str_c("LINE1.full_length.RNAi_2019_paper", "consensus.fasta", sep = ".")))
