### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/05_overlap_with_5pUTR_LINE1s")

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
library(BSgenome.Maur.UCSC.Siomi)
library(systemPipeR)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# LINE1s with complete ORFs path
orfs_path <- file.path(inpath, "L1_full_length_manual_200707.consensus.ORFs_blast_hits.csv")

# LINE1s with long 5pUTRs coordinates
utrs_path <- file.path(inpath, "L1_full_length_manual_200707.manual_complete_5pUTRs.200708.coordinates.csv")

######################################################## READ DATA
# read ORFs table
orfs_tb <- readr::read_csv(orfs_path)

# read 5pUTRs table
utrs_tb <- readr::read_csv(utrs_path)

######################################################## MAIN CODE
# clean ORFs table
orfs_tb %<>% 
  dplyr::select(hit_coordinates, strand, rmsk_id, repName) %>% 
  dplyr::mutate(repName = stringr::str_split(repName, "/") %>% purrr::map(., unique) %>% unlist(.)) %>% 
  tidyr::separate(hit_coordinates, c("seqnames", "start", "end"), sep = " ")

# create GenomicRanges
orfs_gr <- GenomicRanges::GRanges(orfs_tb)
utrs_gr <- GenomicRanges::GRanges(utrs_tb)

# overlap
overlaps <- findOverlaps(utrs_gr, orfs_gr, ignore.strand = T, minoverlap = 5000)

# for each LINE1 with proper 5'UTRs add info about rmskID and repName
utrs_gr_hits <- utrs_gr[queryHits(overlaps)]
orfs_gr_hits <- orfs_gr[subjectHits(overlaps)]
mcols(utrs_gr_hits)$rmsk_id <- mcols(orfs_gr_hits)$rmsk_id
mcols(utrs_gr_hits)$repName <- mcols(orfs_gr_hits)$repName
strand(utrs_gr_hits) <- strand(orfs_gr_hits)


### find ORFs 
# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Maur.UCSC.Siomi, names = utrs_gr_hits)
names(line1_seq) <- mcols(utrs_gr_hits)$rmsk_id

# find ORFs
line1_orfs <- predORF(line1_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = F)

# join to one GRanges object
insertions_orfs_list <- 
  line1_orfs %>% 
  unname(.) %>% 
  do.call(c, .)

# set repeatMasker id
mcols(insertions_orfs_list)$rmsk_id <- as.character(seqnames(insertions_orfs_list))

# split by name and reading frame
insertions_orfs_list <- split(insertions_orfs_list, list(insertions_orfs_list$rmsk_id, insertions_orfs_list$inframe2end))


### reduce ORF coordinates per reading frame and strand in each insertion
# loop
insertions_orfs_reduced <- purrr::map(names(insertions_orfs_list), function(orf_name){
  
  # get one insertion/frame
  insertions_orf <- insertions_orfs_list[[orf_name]]
  
  # get name and reading frame
  rmsk_id <- unique(mcols(insertions_orf)$rmsk_id)
  
  # to GRanges, reduce
  insertions_orf <- 
    insertions_orfs_list[[orf_name]] %>% 
    GenomicRanges::reduce(., ignore.strand = F, min.gapwidth = 0)
  
  # set rmsk_id and reading frames
  mcols(insertions_orf)$rmsk_id <- rmsk_id
  
  # return 
  return(insertions_orf)
  
})

### find two longest ORFs per rmsk_id and save
# get tidy table
insertions_orfs_tb <- 
  insertions_orfs_reduced %>% 
  do.call(c,.) %>% 
  as_tibble(.)

insertions_orfs_tb %<>% 
  dplyr::select(rmsk_id, start, end, strand, width) %>%
  dplyr::arrange(-width) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate_at(vars(dplyr::starts_with("longest_orf")), ~replace(., is.na(.), 0)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(longest_orf_1 >= 1200*3, 
                longest_orf_2 >= 370*3)

# add info about ORF to final table
line1_tb <- 
  utrs_gr_hits %>% 
  as_tibble(.) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::left_join(., insertions_orfs_tb, by = "rmsk_id") %>% 
  dplyr::select(-query_id)

# save as table
readr::write_csv(line1_tb, file.path(outpath, "LINE1.complete_ORFs.long_5pUTR.20200729.csv"))

# # save as SAF
# line1_tb %>% 
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>% 
#   readr::write_delim(., file.path(outpath, "LINE1.complete_ORFs.long_5pUTR.20200729.saf"), delim = "\t")

### save as .bed
# to GRanges
line1_bed <- 
  line1_tb %>% 
  GRanges(.)

# set names
names(line1_bed) <- mcols(line1_bed)$rmsk_id

# save
rtracklayer::export.bed(line1_bed, file.path(outpath, "LINE1.complete_ORFs.long_5pUTR.20200729.bed"))


