### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/rat.rn6/LINE1/substitution_rate/random_sampled_LINE1s")

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

library(BSgenome.Rnorvegicus.UCSC.rn6)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/rat/rn6.Rnor_6.0.GCA_000001895.4"

# raw repeatMasker
rmsk_raw_path <- list.files(genome_dir, "rmsk\\..*\\.raw\\.fa\\.out\\.gz", full.names = T)

######################################################## READ DATA
# read raw repeatMasker
rmsk_raw <- readr::read_table2(file = rmsk_raw_path, skip = 3, 
                               col_names = c("bit_score", "perc_substitution", "perc_deletion", "perc_insertion", 
                                             "seqnames", "query_start", "query_end", "query_left", "strand",
                                             "repName", "repClass_repFamily", "repeat_begin", "repeat_start", "repeat_end", 
                                             "rmsk_id", "tmp"))

######################################################## MAIN CODE
### filter raw repeatMasker table
# get width and percent of alignment of all elements
line1_raw <- 
  rmsk_raw %>% 
  dplyr::select(seqnames, start = query_start, end = query_end, strand, 
                bit_score, perc_substitution, perc_deletion, perc_insertion,
                repName, repClass_repFamily, 
                repeat_begin, repeat_start, repeat_end, 
                rmsk_id) %>% 
  dplyr::filter(repClass_repFamily == "LINE/L1") %>% 
  dplyr::mutate(strand = str_replace(strand, "C", "-")) %>% 
  dplyr::mutate(repeatStart = ifelse(strand == "+", repeat_begin, repeat_end), 
                repeatEnd = repeat_start, 
                repeatLeft = ifelse(str_detect(repeat_begin, "\\("), repeat_begin, repeat_end)) %>% 
  dplyr::mutate(repeatStart = as.numeric(repeatStart),
                repeatEnd = as.numeric(repeatEnd), 
                repeatLeft = str_remove_all(repeatLeft, "\\(|\\)") %>% as.numeric(.)) %>% 
  dplyr::mutate(repeatWidth = repeatEnd - repeatStart + 1, 
                genomeWidth = end - start + 1,
                consensusWidth = repeatStart + repeatWidth + repeatLeft - 1) %>% 
  dplyr::select(seqnames:strand, repName, repeatStart:repeatLeft, repeatWidth, genomeWidth, consensusWidth, rmsk_id) 

# calculate overlap with consensus sequence
line1_filt <- 
  line1_raw %>% 
  dplyr::mutate(perc_align = 100 * round(genomeWidth / consensusWidth, 3))
                                                                                                                                                                                    
### sample 200 LINE1s from each repName
# sample
set.seed(1234)
line1_sample_tb <- 
  line1_filt %>%
  dplyr::group_by(repName) %>% 
  sample_n(ifelse(n() >= 200, 200, n())) %>% 
  dplyr::ungroup(.)

# save table
readr::write_csv(line1_sample_tb, file.path(outpath, "LINE1s.random_200_per_repName.200806.csv"))


### get sequences
# to GRanges
line1_sample_gr <- 
  line1_sample_tb %>% 
  GRanges(.) 

# get sequences 
line1_sample_seq <- Biostrings::getSeq(BSgenome.Rnorvegicus.UCSC.rn6, line1_sample_gr)


### write fasta by subfamilies
# split by subfamily
line1_sample_seq_list <- base::split(line1_sample_seq, mcols(line1_sample_gr)$repName)

# write
purrr::map(names(line1_sample_seq_list), function(rep_name){
  
  # get one LTR subfamily, save
  line1_sample_seq_list[[rep_name]] %T>% 
    Biostrings::writeXStringSet(., file.path(outpath, "fasta_files", str_c("LINE1s.random_200_per_repName", rep_name, "200806.single.fasta", sep = ".")))
  
  # return
  return(rep_name)
  
})



