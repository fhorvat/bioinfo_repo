### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate")

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

# repeatMasker with joined consecutive rmsk ID
rmsk_joined_path <- file.path(inpath, "rmsk.Siomi.20200701.joined_consecutive.fa.out.gz")

# chosen LTR classes path
documentation_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate"
classes_tb_path <- file.path(documentation_path, "chosen_LTR_classes.csv")

######################################################## READ DATA
# read repeatMasker with joined consecutive rmsk ID
rmsk_joined_tb <- readr::read_delim(rmsk_joined_path, delim = "\t")

# read table with chosen LTR classes
classes_tb <- readr::read_csv(classes_tb_path)
  
######################################################## MAIN CODE
### sample 200 LTRs of full length
# get all individual repeatMasker insertions which were not interupted and which are >= 80% of the consensus sequence
ltr_individual_tb <-
  rmsk_joined_tb %>%
  dplyr::right_join(., classes_tb, by = "repName") %>% 
  dplyr::mutate(whole_width = repeatStart + repeatWidth + repeatLeft - 1, 
                genome_width = end - start + 1, 
                perc_align = 100 * round(genome_width / whole_width, 3)) %>% 
  dplyr::filter(perc_align >= 80)
  
# sample LTRs - 200 in each repeat subfamily
set.seed(1234)
ltr_sample_tb <- 
  ltr_individual_tb %>%
  dplyr::group_by(repSubfamily) %>% 
  sample_n(ifelse(n() >= 200, 200, n())) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(repSubfamily = factor(repSubfamily, levels = names(comboClasses))) %>% 
  dplyr::arrange(repSubfamily)

# save table
readr::write_csv(ltr_sample_tb, file.path(outpath, "LTRs.random_200_per_subfamily.80_percent_of_length.csv"))


### get sequences
# to GRanges
ltr_sample_gr <- 
  ltr_sample_tb %>% 
  GRanges(.) 

# get sequences 
ltr_sample_seq <- Biostrings::getSeq(BSgenome.Maur.UCSC.Siomi, ltr_sample_gr)

# save all
Biostrings::writeXStringSet(ltr_sample_seq, file.path(outpath, "LTRs.random_200_per_subfamily.80_percent_of_length.all.fasta"))


### write fasta by subfamilies
# split by subfamily
ltr_sample_seq_list <- base::split(ltr_sample_seq, mcols(ltr_sample_gr)$repSubfamily)

# write
purrr::map(names(ltr_sample_seq_list), function(ltr_name){
  
  # get one LTR subfamily, save
  ltr_sample_seq_list[[ltr_name]] %T>% 
    Biostrings::writeXStringSet(., file.path(outpath, str_c("LTRs.random_200_per_subfamily", ltr_name, "80_percent_of_length.single.fasta", sep = ".")))
  
  # return
  return(ltr_name)
  
})



