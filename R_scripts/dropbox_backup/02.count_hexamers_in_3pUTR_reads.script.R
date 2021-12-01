#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts hexamers in 3pUTRs 
### DATE: Thu Sep 16 01:31:05 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# wideScreen()

######################################################## WORKING DIRECTORY
setwd(".")

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
library(GenomicAlignments)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
bam_path <- args$bam_path
bam_name <- args$bam_name
utr_info_path <- args$utr_info_path
bsgenome <- args$bsgenome

# load BSGenome package
lapply(bsgenome, library, character.only = TRUE)

######################################################## READ DATA
# read coordinates
utr_coords <- readr::read_csv(utr_info_path)

######################################################## MAIN CODE
# get GRanges
utr_gr <- GRanges(utr_coords)

# reduce 3'UTR coordinates
utr_gr_reduced <- reduce(utr_gr, ignore.strand = T)

# overlap with 3' UTRs
overlap <- findOverlaps(utr_gr_reduced, utr_gr, ignore.strand = T)

# set strand
utr_gr_reduced_hits <- 
  utr_gr_reduced[queryHits(overlap)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(hit_id = queryHits(overlap), 
                strand_utr = as.character(strand(utr_gr[subjectHits(overlap)]))) %>% 
  dplyr::group_by(hit_id) %>% 
  dplyr::summarise(strand_utr = str_c(unique(strand_utr), collapse = "|")) %>% 
  dplyr::filter(!str_detect(strand_utr, "\\|"))

# get UTR coordinates which don't overlap both strands
utr_gr_reduced <- utr_gr_reduced[utr_gr_reduced_hits$hit_id]
strand(utr_gr_reduced) <- utr_gr_reduced_hits$strand_utr


### assign reads to 3'UTR
# create empty vector 
kmer_count_total <- rep(0, 4096)

# open connection to bam file in chunks
bamfile <- BamFile(bam_path, yieldSize = 100000)
open(bamfile)

# load chunks of alignments from bam file and classify each alignment
while(length(bam_gr_list <- readGAlignmentsList(bamfile, 
                                                param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = F))))){
  
  # get read coordinates
  bam_gr <- 
    bam_gr_list %>% 
    grglist(.) %>% 
    reduce(., ignore.strand = T) %>% 
    unlist(.)
  
  # overlap
  overlaps <- findOverlaps(bam_gr, utr_gr_reduced, ignore.strand = T)
  
  # intersect with 3'UTRs to get only coordinates which really overlap 3'UTR
  bam_gr_utr <- bam_gr[queryHits(overlaps)] 
  bam_gr_utr <- pintersect(bam_gr_utr, utr_gr_reduced[subjectHits(overlaps)], ignore.strand = T)
  
  # set strand to read
  strand(bam_gr_utr) <- strand(utr_gr_reduced[subjectHits(overlaps)])
  
  # get sequences 
  bam_gr_sequences <- Biostrings::getSeq(get(bsgenome), bam_gr_utr)
  
  # count hexamers
  kmer_count <- 
    Biostrings::oligonucleotideFrequency(bam_gr_sequences, 6, as.prob = F, with.labels = TRUE) %>% 
    colSums(.)
  
  # add to vector
  kmer_count_total <- kmer_count_total + kmer_count
  
}

# create table for plot
kmer_freq_tb <- 
  tibble(kmer = names(kmer_count_total), 
         count = kmer_count_total) %>% 
  dplyr::mutate(freq = 100*(count / sum(count))) %>% 
  dplyr::arrange(-freq) %>% 
  dplyr::mutate(pos = 1:n(),
                pos = factor(pos, levels = 1:n()))

# save table
readr::write_csv(kmer_freq_tb, file.path(outpath, str_c(bam_name, "6mer_counts", "csv", sep = ".")))

# plot as points
kmer_points <- 
  ggplot(kmer_freq_tb, aes(x = pos, y = freq), fill = "black") + 
  geom_point() +
  ylab("") + 
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save 
ggsave(filename = file.path(outpath, str_c(bam_name, "6mer_counts", "png", sep = ".")), 
       plot = kmer_points, 
       width = 10, height = 10)
