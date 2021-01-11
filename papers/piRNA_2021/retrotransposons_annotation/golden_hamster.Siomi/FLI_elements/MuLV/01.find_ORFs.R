### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MuLV/annotation")

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
rmsk_joined_path <- file.path(genome_dir, "rmsk.Siomi.20200701.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

# clean repeatModeler path
rmod_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/RepeatModeler/RepeatMasker"
rmod_clean_path <- file.path(rmod_path, "rmsk.Siomi.20200728.RepeatModeler.clean.fa.out.gz")
rmod_joined_path <- file.path(rmod_path, "rmsk.Siomi.20200728.RepeatModeler.joined_rmsk_id.fa.out.gz")

######################################################## READ DATA
# read joined repeatMasker
rmsk_joined <- 
  readr::read_delim(rmsk_path, delim = "\t") %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

# read clean repeatMasker
rmsk_clean <- 
  readr::read_delim(rmsk_clean_path, delim = "\t") %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

# read clean repeatModeler
rmod_clean <- 
  readr::read_delim(rmod_clean_path, delim = "\t") %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

# read repeatModeler
rmod_joined <- 
  readr::read_delim(rmod_joined_path, delim = "\t") %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

######################################################## MAIN CODE
# filter clean repeatMasker
rmsk_clean_filt <- 
  rmsk_clean %>%
  dplyr::filter(repClass == "LTR") %>% 
  dplyr::filter(str_detect(repName, "MuLV")) %>% 
  dplyr::arrange(-width) %>% 
  tidyr::unite(col = "coordinates", seqnames, start, end, sep = " ")

# # save
# readr::write_csv(rmsk_clean_filt, file.path(outpath, "rmsk.Siomi.20200701.MuLV.coordinates.csv"))

### overlap with repeatModeler coordinates
# create repeatModeler GRanges
rmod_gr <- 
  rmod_joined %>% 
  dplyr::mutate(strand = ifelse(strand %in% c("+", "-"), strand, "*")) %>% 
  GRanges(.)

# create repeatMasker GRanges
rmsk_gr <- 
  rmsk_clean_filt %>%
  tidyr::separate(coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# overlap
overlap <- findOverlaps(rmsk_gr, rmod_gr, ignore.strand = F)

# get hits in repeatModeler 
rmod_hits <- 
  rmod_gr[subjectHits(overlap)] %>% 
  unique(.) %>% 
  as_tibble(.) %>% 
  dplyr::filter(insertion_class == "whole") %>% 
  GRanges(.) %>% 
  unique(.)

### find ORFs 
# extract sequences
rmod_seq <- Biostrings::getSeq(x = BSgenome.Maur.UCSC.Siomi, names = rmod_hits)
names(rmod_seq) <- mcols(rmod_hits)$rmsk_id

# find ORFs
rmod_orfs <- predORF(rmod_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)

# get tidy table
rmod_orfs_tb<- 
  rmod_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(rmsk_id = seqnames, width) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.)

# add info about ORF to final table
rmod_tb <- 
  rmod_hits %>% 
  as_tibble(.) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::left_join(., rmod_orfs_tb, by = "rmsk_id") %>% 
  arrange(desc(longest_orf_1 + longest_orf_2))

