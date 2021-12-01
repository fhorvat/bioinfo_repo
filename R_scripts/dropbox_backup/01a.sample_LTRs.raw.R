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
library(rtracklayer)
library(DESeq2)

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

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

# raw repeatMasker
rmsk_raw_path <- file.path(genome_dir, "rmsk.Siomi.20200701.raw.fa.out.gz")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")

# read raw repeatMasker
rmsk_raw <- readr::read_table2(file = rmsk_raw_path, skip = 3, col_names = F)

######################################################## MAIN CODE
# IAP, MT, ORR, RMER, MER, MER52, RodERV21
comboClasses <- list("MT2"   = c("MT2_Mm"),
                     "MT2A"  = c("MT2A"),
                     "MT2B"  = c("MT2B", "MT2B1", "MT2B2"),
                     "MT2C"  = c("MT2C_Mm"),
                     "MTA"   = c("MTA_Mm"),
                     "MTB"   = c("MTB", "MTB_Mm"),
                     "MTC"   = c("MTC"),
                     "MTD"   = c("MTD"),
                     "MTE"   = c("MTEa", "MTEb"),
                     "MTE2"  = c("MTE2a", "MTE2b"),
                     "ORR1A" = c("ORR1A0", "ORR1A1", "ORR1A2", "ORR1A3", "ORR1A4"),
                     "ORR1B" = c("ORR1B1", "ORR1B2"),
                     "ORR1C" = c("ORR1C1", "ORR1C2"),
                     "ORR1D" = c("ORR1D1", "ORR1D2"),
                     "ORR1E" = c("ORR1E"),
                     "ORR1F" = c("ORR1F"),
                     "ORR1G" = c("ORR1G"), 
                     "IAP1" = c("IAP1-MM_LTR", "IAP1-MM_I"), 
                     "IAPLTR" = c("IAPLTR1_Mm", "IAPLTR1a_Mm", "IAPLTR2_Mm", "IAPLTR2a", "IAPLTR2a2_Mm", "IAPLTR2b", "IAPLTR3", "IAPLTR4", "IAPLTR4_I"), 
                     "IAPEY" = c("IAPEY_LTR", "IAPEY2_LTR", "IAPEY3_LTR", "IAPEY3C_LTR", "IAPEY4_LTR", "IAPEY4_I", "IAPEY5_I", "IAPEY5_LTR"),
                     "RMER1" = c("RMER10A", "RMER10B", "RMER12", "RMER12B", "RMER12C", "RMER13A", "RMER13A1", "RMER13A2", "RMER13B", "RMER15", "RMER16", 
                                 "RMER16A2", "RMER16A3", "RMER16B", "RMER16B2", "RMER16B3", "RMER16C", "RMER16_Mm", "RMER17A", "RMER17A2", "RMER17B", 
                                 "RMER17B2", "RMER17C", "RMER17C2", "RMER17D", "RMER17D2", "RMER19A", "RMER19B", "RMER19B2", "RMER19C"), 
                     "RMER2x" = c("RMER2", "RMER20A", "RMER20B", "RMER20C_Mm", "RMER21A", "RMER21B"), 
                     "RMER3x" = c("RMER3D1", "RMER3D2", "RMER3D3", "RMER3D4"), 
                     "RMER4x" = c("RMER4A", "RMER4B"),
                     "RMER5x" = c("RMER5"), 
                     "RMER6x" = c("RMER6A", "RMER6B", "RMER6BA", "RMER6C", "RMER6D"), 
                     "MER1x" = c("MER101B", "MER110", "MER110A", "MER129"), 
                     "MER2x" = c("MER21B", "MER21C"), 
                     "MER3x" = c("MER31A", "MER31B", "MER34", "MER34A", "MER34A1", "MER34B", "MER34C", "MER34C_v"), 
                     "MER5x" = c("MER54A", "MER54B", "MER57C2", "MER57D", "MER57E1", "MER57E2"), 
                     "MER6x" = c("MER65D", "MER67A", "MER67B", "MER67C", "MER67D", "MER68", "MER68B", "MER68C"), 
                     "MER7x" = c("MER70A", "MER70B", "MER70C", "MER73", "MER74A", "MER74B", "MER74C", "MER76", "MER77", "MER77B"), 
                     "MER8x" = c("MER88", "MER89"), 
                     "MER9x" = c("MER90", "MER90a", "MER92A", "MER92B", "MER92C", "MER92D", "MER95"))

# list to tibble
classes_tb <- purrr::map(names(comboClasses), function(ltr_name){
  
  # get one LTR type in tibble
  comboClasses[[ltr_name]] %>% 
    tibble(repName = ., repSubfamily = ltr_name)
  
}) %>%
  dplyr::bind_rows(.)

# sample LTRs - 200 in each repeat subfamily
set.seed(1234)
ltr_sample_tb <- 
  rmsk_clean %>%
  dplyr::right_join(., classes_tb, by = "repName") %>% 
  dplyr::group_by(repSubfamily) %>% 
  sample_n(ifelse(n() >= 200, 200, n())) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(repSubfamily = factor(repSubfamily, levels = names(comboClasses))) %>% 
  dplyr::arrange(repSubfamily)

# save table
readr::write_csv(ltr_sample_tb, file.path(outpath, "LTRs.random_200_per_subfamily.csv"))


### get sequences
# to GRanges
ltr_sample_gr <- 
  ltr_sample_tb %>% 
  GRanges(.) 

# get sequences 
ltr_sample_seq <- Biostrings::getSeq(BSgenome.Maur.UCSC.Siomi, ltr_sample_gr)

# save all
Biostrings::writeXStringSet(ltr_sample_seq, file.path(outpath, "LTRs.random_200_per_subfamily.all.fasta"))


### write fasta by subfamilies
# split by subfamily
ltr_sample_seq_list <- base::split(ltr_sample_seq, mcols(ltr_sample_gr)$repSubfamily)

# write
purrr::map(names(ltr_sample_seq_list), function(ltr_name){
  
  # get one LTR subfamily, save
  ltr_sample_seq_list[[ltr_name]] %T>% 
    Biostrings::writeXStringSet(., file.path(outpath, str_c("LTRs.random_200_per_subfamily", ltr_name, "single.fasta", sep = ".")))
  
  # return
  return(ltr_name)
  
})



