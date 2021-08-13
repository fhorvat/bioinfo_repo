### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/annotation/01_rmsk_annotation")

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

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

######################################################## MAIN CODE
# get widths
rmsk_tb %<>%    
  dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "+", strand)) %>% 
  GRanges(.) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(rmsk_id = as.character(rmsk_id))

# filter LINE1s
line1_tb <- 
  rmsk_tb %>% 
  dplyr::filter(repFamily == "L1", 
                width > 4000)

# to GRanges
line1_gr <- GRanges(line1_tb) 

# extract sequences
line1_seq <- Biostrings::getSeq(x = BSgenome.Maur.UCSC.Siomi, names = line1_gr)
names(line1_seq) <- line1_gr$rmsk_id

# find sequences containing Ns
n_sequences <- 
  vmatchPattern(pattern = "N", subject = line1_seq) %>% 
  unlist(.) %>% 
  names(.) %>% 
  unique(.)

# remove sequences containg N
line1_seq <- line1_seq[!names(line1_seq) %in% n_sequences]

# find ORFs
# line1_orfs <- predORF(line1_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = F)
# saveRDS(line1_orfs, file = file.path(outpath, "LINE1.4000nt_plus.ORFs.grl.RDS"))
line1_orfs <- readRDS(file = file.path(outpath, "LINE1.4000nt_plus.ORFs.grl.RDS"))
  
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
  dplyr::mutate_at(vars(dplyr::starts_with("longest_orf")), ~replace(., is.na(.), 0))

# add info about ORFs to original table, save
line1_tb %<>% 
  dplyr::left_join(., insertions_orfs_tb, by = "rmsk_id") %>% 
  dplyr::filter(rmsk_id %in% names(line1_seq)) %T>%
  readr::write_csv(., file.path(outpath, "LINE1.4000nt_plus.ORFs.csv"))

