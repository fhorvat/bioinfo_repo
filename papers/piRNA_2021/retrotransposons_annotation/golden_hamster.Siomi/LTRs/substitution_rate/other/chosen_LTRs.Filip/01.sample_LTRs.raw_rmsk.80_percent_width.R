### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate/random_sampled_LTRs")

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

library(BSgenome.Maur.UCSC.Siomi)
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
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# raw repeatMasker
rmsk_raw_path <- file.path(genome_dir, "rmsk.Siomi.20200701.raw.fa.out.gz")

# chosen LTR classes path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate"
classes_tb_path <- file.path(documentation_path, "chosen_LTR_classes.csv")

######################################################## READ DATA
# read raw repeatMasker
rmsk_raw <- readr::read_table2(file = rmsk_raw_path, skip = 3, 
                               col_names = c("bit_score", "perc_substitution", "perc_deletion", "perc_insertion", 
                                             "seqnames", "query_start", "query_end", "query_left", "strand",
                                             "repName", "repClass_repFamily", "repeat_begin", "repeat_start", "repeat_end", 
                                             "rmsk_id", "tmp"))

# read table with chosen LTR classes
classes_tb <- readr::read_csv(classes_tb_path)

######################################################## MAIN CODE
### filter raw repeatMasker table
# get width and percent of alignment of all LTRs
ltrs_raw <- 
  rmsk_raw %>% 
  dplyr::select(seqnames, start = query_start, end = query_end, strand, 
                bit_score, perc_substitution, perc_deletion, perc_insertion,
                repName, repClass_repFamily, 
                repeat_begin, repeat_start, repeat_end, 
                rmsk_id) %>% 
  dplyr::right_join(., classes_tb, by = "repName") %>% 
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
  dplyr::select(seqnames:strand, repName, repSubfamily, repeatStart:repeatLeft, repeatWidth, genomeWidth, consensusWidth, rmsk_id) 

# get LTRs with width >= 80% of consensus sequence
ltrs_filt <- 
  ltrs_raw %>% 
  dplyr::mutate(perc_align = 100 * round(genomeWidth / consensusWidth, 3)) %>% 
  dplyr::filter(perc_align >= 80)


### sample 200 LTRs of full length
# sample LTRs - 200 in each repeat subfamily
set.seed(1234)
ltr_sample_tb <- 
  ltrs_filt %>%
  dplyr::group_by(repSubfamily) %>% 
  sample_n(ifelse(n() >= 200, 200, n())) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(repSubfamily = factor(repSubfamily, levels = unique(classes_tb$repSubfamily))) %>% 
  dplyr::arrange(repSubfamily)

# save table
readr::write_csv(ltr_sample_tb, file.path(outpath, "LTRs.random_200_per_subfamily.80_percent_of_length.raw_rmsk.csv"))


### get sequences
# to GRanges
ltr_sample_gr <- 
  ltr_sample_tb %>% 
  GRanges(.) 

# get sequences 
ltr_sample_seq <- Biostrings::getSeq(BSgenome.Maur.UCSC.Siomi, ltr_sample_gr)

# # save all
# Biostrings::writeXStringSet(ltr_sample_seq, file.path(outpath, "LTRs.random_200_per_subfamily.80_percent_of_length.all.fasta"))


### write fasta by subfamilies
# split by subfamily
ltr_sample_seq_list <- base::split(ltr_sample_seq, mcols(ltr_sample_gr)$repSubfamily)

# write
purrr::map(names(ltr_sample_seq_list), function(ltr_name){
  
  # get one LTR subfamily, save
  ltr_sample_seq_list[[ltr_name]] %T>% 
    Biostrings::writeXStringSet(., file.path(outpath, str_c("LTRs.random_200_per_subfamily", ltr_name, "80_percent_of_length.raw_rmsk.single.fasta", sep = ".")))
  
  # return
  return(ltr_name)
  
})



