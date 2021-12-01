### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/MYSERV.repeatMasker/annotation")

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
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# annotation table path
annotation_path <- file.path(inpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.csv")

######################################################## READ DATA
# read annotation
annotation_tb <- readr::read_csv(annotation_path)

######################################################## MAIN CODE
### get insertion sequences
# get ranges
annotation_gr <- 
  annotation_tb %>% 
  tidyr::separate(col = coordinates, into = c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
  GRanges(.)

# get sequences
insertions_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, annotation_gr)
names(insertions_seq) <- mcols(annotation_gr)$coordinates

# find sequences containing Ns
n_sequences <- 
  vmatchPattern(pattern = "N", subject = insertions_seq) %>% 
  unlist(.) %>% 
  names(.) %>% 
  unique(.)

# remove sequences containg N
insertions_seq <- insertions_seq[!names(insertions_seq) %in% n_sequences]


### find ORFs
# find all ORFs on both strands
insertions_orfs <- predORF(insertions_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = F)

# join to one GRanges object
insertions_orfs_list <- 
  insertions_orfs %>% 
  unname(.) %>% 
  do.call(c, .)

# set coordinates
mcols(insertions_orfs_list)$coordinates <- as.character(seqnames(insertions_orfs_list))

# split by name and reading frame
insertions_orfs_list <- split(insertions_orfs_list, list(insertions_orfs_list$coordinates, insertions_orfs_list$inframe2end))


### reduce ORF coordinates per reading frame and strand in each insertion
# loop
insertions_orfs_reduced <- purrr::map(names(insertions_orfs_list), function(orf_name){
  
  # get one insertion/frame
  insertions_orf <- insertions_orfs_list[[orf_name]]
  
  # get name and reading frame
  coordinates <- unique(mcols(insertions_orf)$coordinates)
  
  # to GRanges, reduce
  insertions_orf <- 
    insertions_orfs_list[[orf_name]] %>% 
    GenomicRanges::reduce(., ignore.strand = F, min.gapwidth = 0)
  
  # set coordinates and reading frames
  mcols(insertions_orf)$coordinates <- coordinates
  
  # return 
  return(insertions_orf)
  
})


### find two longest ORFs per rmsk_id and save
# get tidy table
insertions_orfs_tb <- 
  insertions_orfs_reduced %>% 
  do.call(c,.) %>% 
  as_tibble(.) %>% 
  dplyr::select(coordinates, start, end, strand, width) %>%
  dplyr::arrange(-width) %>% 
  dplyr::group_by(coordinates) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate_at(vars(dplyr::starts_with("longest_orf")), ~replace(., is.na(.), 0))

# add info about ORFs to original table, filter
annotation_tb_full <- 
  annotation_tb %>% 
  dplyr::left_join(., insertions_orfs_tb, by = "coordinates") %>% 
  dplyr::arrange(desc(longest_orf_1 + longest_orf_2))

# save 
readr::write_csv(annotation_tb_full, file.path(outpath, "MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.csv"))


