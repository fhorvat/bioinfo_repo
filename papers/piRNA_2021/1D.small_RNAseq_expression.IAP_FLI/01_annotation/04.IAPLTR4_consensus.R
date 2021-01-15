### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/annotation/IAPLTR4_consensus")

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
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(stringdist)
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

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# insertions table path
insertions_path <- file.path(inpath, "..", "IAP.potentially_young.ordered_by_ORFs.20201031.top_110.csv")

# LTRs table path
ltrs_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/potentially_young_LTRs"
ltrs_path <- file.path(ltrs_path, "potentially_young_LTRs.edit_distance.10.ORFs.csv")

######################################################## READ DATA
# read insertion table
insertions_tb <- readr::read_csv(insertions_path)

# read LTRs table
ltrs_tb <- readr::read_csv(ltrs_path)

######################################################## MAIN CODE
# clean LTRs table
ltrs_tb %<>% 
  dplyr::select(rmsk_id, repName, LTRs_align_score = align_score, LTRs_edit_distance = edit_distance, 
                coordinates_5pLTR, coordinates_3pLTR, sequence_5pLTR, sequence_3pLTR)

# filter 
iap_tb <- 
  insertions_tb %>% 
  dplyr::filter(ltr_name == "IAPLTR4_I") %>% 
  dplyr::filter(longest_orf_1 == 2370, longest_orf_2 == 1821) %>% 
  dplyr::left_join(., ltrs_tb, by = "rmsk_id")

# save table
iap_tb %>% 
  tidyr::unite("coordinates", seqnames, start, end, sep = " ") %>% 
  dplyr::select(coordinates, strand, ltr_name, insertion_name = repName, everything()) %>% 
  dplyr::rename(insertion_width = width) %>% 
  readr::write_csv(., file = file.path(outpath, str_c("IAPLTR4_I", "consistent_ORFs", "csv", sep = ".")))

# to GRanges
iap_gr <-
  iap_tb %>%
  GRanges(.)
names(iap_gr) <- str_c(mcols(iap_gr)$ltr_name, mcols(iap_gr)$rmsk_id, sep = ".")
names(iap_gr) <- str_replace_all(names(iap_gr), "-", "_")

# get sequences
iap_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, iap_gr)

# # save as fasta
Biostrings::writeXStringSet(x = iap_seq, filepath = file.path(outpath, str_c("IAPLTR4_I", "consistent_ORFs", "fasta", sep = ".")))


### get consensus by repeat names
# split by repeat names
iap_seq_list <- split(iap_seq, iap_gr$ltr_name)

## MSA and get consesus
iap_seq_consensus <- purrr::map(names(iap_seq_list), function(seq_name){
  
  # do the MSA
  seq_msa <- msa::msa(iap_seq_list[[seq_name]], method = "ClustalOmega")
  # seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq_list[[seq_name]], iterations = 100, refinements = 100, processors = 1, verbose = T)
  
  # write as fasta
  seq_msa@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, str_c("IAPLTR4_I", "consistent_ORFs", "msa", "fasta", sep = ".")))
  
  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)
  
  # return
  return(seq_consensus)
  
})

# join together
iap_seq_consensus <- do.call(c, iap_seq_consensus)
names(iap_seq_consensus) <- names(iap_seq_list)

# save the consensus
Biostrings::writeXStringSet(iap_seq_consensus, file = file.path(outpath, str_c("IAPLTR4_I", "consistent_ORFs", "consensus", "fasta", sep = ".")))


