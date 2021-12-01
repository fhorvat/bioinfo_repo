### INFO: 
### DATE: Tue Jul 09 19:28:50 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_locus/flanking_genes.mRNA/results/clustalW_3prime")

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

library(rtracklayer)
library(Biostrings)
library(msa)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# lnc1 3' end fasta path
lnc1_3prime_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/rodents_evolution/lnc1_3prime_end/mm10.lnc1_3prime.fa"

# list lnc1 locus fasta files
lnc1_locus_path <- list.files(inpath, ".*\\.inter_exons.fa", full.names = T)

######################################################## READ DATA
# read lnc1 3' end fasta
lnc1_3prime_seq <- 
  Biostrings::readDNAStringSet(lnc1_3prime_path) %>% 
  Biostrings::reverseComplement(.)
  
# read lnc1 locus fasta
lnc1_locus_seq <- Biostrings::readDNAStringSet(lnc1_locus_path)

######################################################## MAIN CODE
# add species to names
names(lnc1_locus_seq) <- str_c(basename(lnc1_locus_path) %>% str_remove("\\.inter_exons\\.fa"), 
                               names(lnc1_locus_seq), sep = ".")

# ### MSA
# # set one lnc1 locus as query
# lnc1_query_name <- "AcoCah"
# lnc1_query <- lnc1_locus_seq[str_detect(names(lnc1_locus_seq), lnc1_query_name)]
# lnc1_query_subseq <- lnc1_query
# # lnc1_query_subseq <- subseq(lnc1_query, 14462, 14734)
# lnc1_query_subseq <- subseq(lnc1_query, 14442, 14754)
# 
# # get sequences of lnc1 locus in genome
# lnc1_pairwise_seq <- c(lnc1_3prime_seq, lnc1_query_subseq)
# 
# # do MSA using ClustalW
# lnc1_msa <- msa::msaClustalW(inputSeqs = lnc1_pairwise_seq)
# 
# # write as fasta
# lnc1_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("lnc1_locus", lnc1_query_name, "lnc1_3prime_mm10", "ClustalW.msa.fasta", sep = ".")))
# 
# # write as .pdf
# msaPrettyPrint(lnc1_msa,
#                output = "pdf",
#                logoColors = "chemical",
#                paperWidth = 32, paperHeight = 5,
#                file = file.path(outpath, str_c("lnc1_locus", lnc1_query_name, "lnc1_3prime_mm10", "ClustalW.msa.pdf", sep = ".")),
#                askForOverwrite = F)


# ### MSA
# # set one lnc1 locus as query
# query_genome <- "MerUng"
# lnc1_query_name <- names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), query_genome)]
# 
# lnc1_query <-
#   lnc1_locus_seq[lnc1_query_name] %>%
#   as.character(.) %>%
#   str_remove_all(., "N") %>%
#   Biostrings::DNAStringSet(.)
# names(lnc1_query) <- lnc1_query_name
# # lnc1_query_subseq <- subseq(lnc1_query, 16850, 17115)
# lnc1_query_subseq <- subseq(lnc1_query, 16830, 17135)
# 
# # get sequences of lnc1 locus in genome
# lnc1_pairwise_seq <- c(lnc1_3prime_seq, lnc1_query_subseq)
# 
# # do MSA using ClustalW
# lnc1_msa <- msa::msaClustalW(inputSeqs = lnc1_pairwise_seq)
# 
# # write as fasta
# lnc1_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("lnc1_locus", query_genome, "lnc1_3prime_mm10", "ClustalW.msa.fasta", sep = ".")))
# 
# # write as .pdf
# msaPrettyPrint(lnc1_msa,
#                output = "pdf",
#                logoColors = "chemical",
#                paperWidth = 32, paperHeight = 5,
#                file = file.path(outpath, str_c("lnc1_locus", query_genome, "lnc1_3prime_mm10", "ClustalW.msa.pdf", sep = ".")),
#                askForOverwrite = F)

# ### MSA
# # set one lnc1 locus as query
# query_genome <- "rn6"
# lnc1_query_name <- names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), query_genome)]
# 
# lnc1_query <-
#   lnc1_locus_seq[lnc1_query_name] %>%
#   as.character(.) %>%
#   str_remove_all(., "N") %>%
#   Biostrings::DNAStringSet(.)
# names(lnc1_query) <- lnc1_query_name
# # lnc1_query_subseq <- subseq(lnc1_query, 6622, 6860)
# lnc1_query_subseq <- subseq(lnc1_query, 6602, 6880)
# 
# # get sequences of lnc1 locus in genome
# lnc1_pairwise_seq <- c(lnc1_3prime_seq, lnc1_query_subseq)
# 
# # do MSA using ClustalW
# lnc1_msa <- msa::msaClustalW(inputSeqs = lnc1_pairwise_seq)
# 
# # write as fasta
# lnc1_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("lnc1_locus", query_genome, "lnc1_3prime_mm10", "ClustalW.msa.fasta", sep = ".")))
# 
# # write as .pdf
# msaPrettyPrint(lnc1_msa,
#                output = "pdf",
#                logoColors = "chemical",
#                paperWidth = 32, paperHeight = 5,
#                file = file.path(outpath, str_c("lnc1_locus", query_genome, "lnc1_3prime_mm10", "ClustalW.msa.pdf", sep = ".")),
#                askForOverwrite = F)


# ### MSA
# # set one lnc1 locus as query
# query_genome <- "ApoSyl"
# lnc1_query_name <- names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), query_genome)]
# 
# lnc1_query <-
#   lnc1_locus_seq[lnc1_query_name] %>%
#   as.character(.) %>%
#   str_remove_all(., "N") %>%
#   Biostrings::DNAStringSet(.)
# names(lnc1_query) <- lnc1_query_name
# # lnc1_query_subseq <- subseq(lnc1_query, 4346, 4612)
# lnc1_query_subseq <- subseq(lnc1_query, 4326, 4632)
# 
# # get sequences of lnc1 locus in genome
# lnc1_pairwise_seq <- c(lnc1_3prime_seq, lnc1_query_subseq)
# 
# # do MSA using ClustalW
# lnc1_msa <- msa::msaClustalW(inputSeqs = lnc1_pairwise_seq)
# 
# # write as fasta
# lnc1_msa@unmasked %>%
#   Biostrings::writeXStringSet(., file = file.path(outpath, str_c("lnc1_locus", query_genome, "lnc1_3prime_mm10", "ClustalW.msa.fasta", sep = ".")))
# 
# # write as .pdf
# msaPrettyPrint(lnc1_msa,
#                output = "pdf",
#                logoColors = "chemical",
#                paperWidth = 32, paperHeight = 5,
#                file = file.path(outpath, str_c("lnc1_locus", query_genome, "lnc1_3prime_mm10", "ClustalW.msa.pdf", sep = ".")),
#                askForOverwrite = F)


### MSA together
# set species order
species_order <- c("mm10", "AcoCah", "MerUng", "rn6", "ApoSyl")

# AcoCah
lnc1_AcoCah <- 
  lnc1_locus_seq[names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "AcoCah")]] %>% 
  subseq(., 14442, 14754)
names(lnc1_AcoCah) <- 
  names(lnc1_AcoCah)[str_detect(names(lnc1_locus_seq), "AcoCah")] %>% 
  str_replace(., ":.*", "14442-14754")

# MerUng
lnc1_MerUng <- 
  lnc1_locus_seq[names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "MerUng")]] %>% 
  as.character(.) %>% 
  str_remove_all(., "N") %>% 
  Biostrings::DNAStringSet(.) %>% 
  subseq(., 16830, 17135)
names(lnc1_MerUng) <- 
  names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "MerUng")] %>% 
  str_replace(., ":.*", "16830-17135")

# rn6
lnc1_rn6 <- 
  lnc1_locus_seq[names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "rn6")]] %>% 
  as.character(.) %>% 
  str_remove_all(., "N") %>% 
  Biostrings::DNAStringSet(.) %>% 
  subseq(., 6602, 6880)
names(lnc1_rn6) <- 
  names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "rn6")] %>% 
  str_replace(., ":.*", "6602-6880")

# ApoSyl
lnc1_ApoSyl <- 
  lnc1_locus_seq[names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "ApoSyl")]] %>% 
  as.character(.) %>% 
  str_remove_all(., "N") %>% 
  Biostrings::DNAStringSet(.) %>% 
  subseq(., 4326, 4632)
names(lnc1_ApoSyl) <- 
  names(lnc1_locus_seq)[str_detect(names(lnc1_locus_seq), "ApoSyl")] %>% 
  str_replace(., ":.*", "4326-4632")



# get sequences of lnc1 locus in genome
lnc1_pairwise_seq <- c(lnc1_3prime_seq, lnc1_AcoCah, lnc1_MerUng, lnc1_rn6, lnc1_ApoSyl)

# do MSA using ClustalW
lnc1_msa <- msa::msaClustalW(inputSeqs = lnc1_pairwise_seq)

# write as fasta
lnc1_msa@unmasked %>% 
  .[names(.) %>% str_remove_all(., "_lnc1:chr7:46967711-46967960|\\..*") %>% match(., species_order) %>% order(.)] %>% 
  Biostrings::writeXStringSet(., file = file.path(outpath, str_c("lnc1_locus", "all_aligned", "ClustalW.msa.fasta", sep = ".")))

# write as .pdf
msaPrettyPrint(lnc1_msa,
               subset = lnc1_msa@unmasked %>% names(.) %>% str_remove_all(., "_lnc1:chr7:46967711-46967960|\\..*") %>% match(., species_order) %>% order(.), 
               output = "pdf", 
               logoColors = "chemical",
               paperWidth = 32, paperHeight = 5,
               file = file.path(outpath, str_c("lnc1_locus", "all_aligned", "ClustalW.msa.pdf", sep = ".")), 
               askForOverwrite = F)

