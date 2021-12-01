### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation")

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

library(BSgenome.Mmusculus.UCSC.mm10)
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
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz")

# reduced exons path
exons_path <- file.path(genome_dir, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.reducedExons.RDS")

# Zoe's list path
zoe_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/LINE1_expression/Documentation/L1s_nested_ours_20180516.ZJM.tidy.csv"

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read Zoe's list
zoe_tb <- readr::read_csv(zoe_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# create GRanges from RepeatMasker
rmsk_gr <- rmsk_tb %>% dplyr::mutate(strand = "*") %>% GRanges(.)

# create GRanges from Zoe
zoe_gr <- GRanges(zoe_tb)


### find two longest ORFs in each annotated LINE1
# extract sequences
zoe_seq <- Biostrings::getSeq(x = BSgenome.Mmusculus.UCSC.mm10, names = zoe_gr)
names(zoe_seq) <- zoe_gr$id

# find ORFs
zoe_orfs <- predORF(zoe_seq, n = "all", type = "grl", mode = "orf", strand = "both", longest_disjoint = T)

# get tidy table
zoe_orfs_tb<- 
  zoe_orfs %>% 
  unlist(.) %>% 
  as_tibble(.) %>% 
  dplyr::select(id = seqnames, width) %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(longest_orf_1 = sort(width, decreasing = T)[1], 
                   longest_orf_2 = sort(width, decreasing = T)[2]) %>% 
  dplyr::ungroup(.) 




# find overlaps and percentage of overlap of each hit
line1_zoe_foverlaps <- findOverlaps(zoe_gr, line1_gr, ignore.strand = T)
overlaps <- pintersect(zoe_gr[queryHits(line1_zoe_foverlaps)], line1_gr[subjectHits(line1_zoe_foverlaps)])
percentOverlap <- width(overlaps) / width(zoe_gr[queryHits(line1_zoe_foverlaps)])

# extract ID's and genes
line1_zoe_tb <- 
  tibble(id = zoe_gr[queryHits(line1_zoe_foverlaps)]$id, 
         rmsk_id = line1_gr[subjectHits(line1_zoe_foverlaps)]$rmsk_id, 
         percent_overlap = percentOverlap)
# dplyr::group_by(id) %>% 
# dplyr::summarise(rmsk_id = str_c(rmsk_id, collapse = ", "))

line1_zoe_tb %>% dplyr::filter(id %in% id[duplicated(id)])

# add rmks info to Zoe's table
zoe_tb_rmsk <- 
  line1_zoe_tb %>% 
  dplyr::left_join(., line1_zoe_tb, by = "id") 


# check if there is any in Zoe's list which doesn't overlap any annotated
# all but one are shorter than 4000, one which is longer is not annotated as one but as two insertions in repeatMasker
zoe_not_present <- 
  zoe_tb %>%
  dplyr::filter(!(id %in% line1_zoe_tb$id)) %>% 
  dplyr::mutate(width = end - start) %>% 
  dplyr::filter(width > 4000)


### filter
# filter based on  ORF length, overlap with annotated exons and integration type
line1_clean <- 
  line1_tb %>% 
  dplyr::filter(longest_orf_1 >= 3800, longest_orf_2 >= 1110) %>% 
  dplyr::filter(is.na(ensembl_id_overlap)) %>% 
  dplyr::filter(insertion_class %in% (c("whole", "within")))


### check Zoe's guys
## prepare
# GRanges
line1_clean_gr <- GRanges(line1_clean)

# find overlaps
zoe_line1_foverlaps <- findOverlaps(zoe_gr, line1_clean_gr, ignore.strand = T)

# extract ID's and genes
zoe_line_tb <- 
  tibble(id = zoe_gr[queryHits(zoe_line1_foverlaps)]$id, 
         rmsk_id = as.character(line1_clean_gr[subjectHits(zoe_line1_foverlaps)]$rmsk_id)) %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(rmsk_id = str_c(rmsk_id, collapse = ", "))

# add to Zoe table
zoe_tb_clean <- 
  zoe_tb %>% 
  dplyr::left_join(., zoe_line_tb, by = "id")


## check in original table
# get Zoe's which don't overlap mine
zoe_missing_id <- zoe_tb_clean %>% dplyr::filter(is.na(rmsk_id)) %$% id

# see what they overlap in original table
zoe_missing_tb <- 
  line1_tb %>% 
  dplyr::filter(zoe_id_overlap %in% zoe_missing_id) %>% 
  dplyr::mutate(orf1_short = (longest_orf_1 < 3800), 
                orf2_short = (longest_orf_2 < 1100)) %>% 
  tidyr::unite(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::select(zoe_id_overlap, coordinates, width, rmsk_id, insertion_class, orf1_short, orf2_short, ensembl_id_overlap) %T>%
  readr::write_csv(., path = file.path(outpath, "missing_zoe.csv"))





## check completely missing ones - ones that don't overlap any LINE1 > 5100
# get ones which don't overlap anything
zoe_completely_missing_id <- zoe_missing_id[!(zoe_missing_id %in% (line1_tb$zoe_id_overlap %>% str_split(., ", ") %>% unlist(.)))]

# filter table
zoe_completely_missing_tb <- 
  zoe_tb %>% 
  dplyr::filter(id %in% zoe_completely_missing_id) %>% 
  GRanges(.) %>% 
  as_tibble(.) %>% 
  dplyr::mutate(is_short = (width < 5100)) %>% 
  dplyr::filter(!is_short)

#### most of them are too short #####

### check remaining 
# GRanges
zoe_completely_missing_gr <- GRanges(zoe_completely_missing_tb)
rmsk_gr <- rmsk_tb %>% dplyr::mutate(strand = "*") %>% GRanges(.)

# overlap with full repeatMasker
zoe_completely_missing_foverlaps <- findOverlaps(zoe_completely_missing_gr, rmsk_gr, ignore.strand = T)

# extract
zoe_completely_missing_rmsk <- rmsk_gr[subjectHits(zoe_completely_missing_foverlaps)]

#### two remaining are classified as more than 1 insertion in repeatMasker so I don't pick them up #### 