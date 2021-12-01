### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate/expression")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# raw repeatMasker
rmsk_raw_path <- file.path(genome_dir, "rmsk.Siomi.20200701.raw.fa.out.gz")

# chosen LTR classes path
classes_tb_path <- file.path(inpath, "..", "all_LTR_classes 200730.xlsx")

######################################################## READ DATA
# read raw repeatMasker
rmsk_raw <- readr::read_table2(file = rmsk_raw_path, skip = 3, 
                               col_names = c("bit_score", "perc_substitution", "perc_deletion", "perc_insertion", 
                                             "seqnames", "query_start", "query_end", "query_left", "strand",
                                             "repName", "repClass_repFamily", "repeat_begin", "repeat_start", "repeat_end", 
                                             "rmsk_id", "tmp"))

# read table with chosen LTR classes
classes_tb <- 
  openxlsx::read.xlsx(classes_tb_path) %>% 
  as_tibble(.) %>% 
  dplyr::select(repName, repSubfamily = category_I)

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
  dplyr::mutate(strand = str_replace(strand, "C", "-")) %>% 
  dplyr::filter(repName %in% classes_tb$repName) %>% 
  dplyr::left_join(., classes_tb, by = "repName")

# format as table
ltrs_coords <- 
  ltrs_raw %>% 
  dplyr::mutate(gene_id = str_c(seqnames, start, end, rmsk_id, sep = ".")) %>% 
  tidyr::unite(col = "coords", seqnames, start, end, sep = " ") %>% 
  dplyr::mutate(repFamily = str_remove(repClass_repFamily, "$LTR/")) %>% 
  dplyr::select(coords, strand, repName, repFamily, repSubfamily, rmsk_id, gene_id)

# save
readr::write_csv(ltrs_coords, file.path(outpath, "rmsk.Siomi.20200701.raw.LTRs.geneInfo.csv"))

# save as SAF
ltrs_coords %>% 
  tidyr::separate(coords, into = c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::select(GeneID = gene_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %T>% 
  readr::write_delim(., file.path(outpath, "rmsk.Siomi.20200701.raw.LTRs.saf"), delim = "\t")


### save as .bed
# to GRanges
ltrs_bed <- 
  ltrs_coords %>% 
  tidyr::separate(coords, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)

# set names
names(ltrs_bed) <- mcols(ltrs_bed)$gene_id

# save
rtracklayer::export.bed(ltrs_bed, file.path(outpath, "rmsk.Siomi.20200701.raw.LTRs.bed"))

