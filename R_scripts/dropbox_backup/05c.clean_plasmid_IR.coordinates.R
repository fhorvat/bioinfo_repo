### INFO: 
### DATE: Mon Oct 01 16:06:52 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_mRNAs.plasmid_coordinates")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# .bed path
bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_plasmid.plasmid_coordinates/ir.plasmid_coordinates.bed"

######################################################## READ DATA
# read bed
ir_gr <- rtracklayer::import(bed_path)

######################################################## MAIN CODE
# clean some stuff
ir_clean <- 
  ir_gr %>% 
  as.data.frame(.) %>%
  as.tibble(.) %>% 
  mutate_if(is.factor, as.character) %>% 
  dplyr::select(-width, -score) %>% 
  # dplyr::filter(!((seqnames == "pCag-EGFP_Lin28IR") & (name == "pCag-EGFP_Lin28IR:3266-3883")),
  #               !((seqnames == "pU6_Lin28IR") & (name == "pCag-EGFP_Lin28IR:3266-3883"))) %>%
  dplyr::mutate(arm = str_extract(seqnames, "Lin28|Mos|mos|Elav2|Elavl2") %>% 
                  str_replace(., "mos", "Mos") %>% 
                  str_replace(., "Elav2", "Elavl2"), 
                arm = ifelse(strand == "+", str_c(arm, ".sense_arm"), str_c(arm, ".antisense_arm"))) %>% 
  dplyr::arrange(seqnames) %T>%
  readr::write_csv(., path = file.path(outpath, "IR_mRNA.plasmid_coordinates.clean.csv"))

# ### read plasmid sequences and check 
# # plasmid paths
# plasmid_path <- 
#   list.files("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/plasmid_sequences", 
#              pattern = "^p.*fasta$", full.names = T) 
# 
# # read plasmids
# plasmid_seq <- Biostrings::readDNAStringSet(plasmid_path)
# 
# # subset plasmids
# plasmid_name <- "pCag-EGFP_Lin28IR" #2854-3167, 3266-3579
# plasmid_subset <- plasmid_seq[plasmid_name]
# 
# # get plus and minus
# int <- 309
# plasmid_plus <- subseq(plasmid_subset, 3266, 3581 + int)
# plasmid_minus <- reverseComplement(subseq(plasmid_subset, 2852 - int, 3167))
# c(plasmid_plus, plasmid_minus)
