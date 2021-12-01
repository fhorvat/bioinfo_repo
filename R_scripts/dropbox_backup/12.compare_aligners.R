### INFO: 
### DATE: Sun Sep 30 00:51:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/IR_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# STAR mapped path
star_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR.only_plasmids"
star_sample_path <- file.path(star_path, "s_pCag_EGFP_Lin28IR_r1.SE.mis_0.bam")

# shrimp mapped path
shrimp_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/shrimp.only_plasmids"
shrimp_sample_path <- file.path(shrimp_path, "s_pCag_EGFP_Lin28IR_r1.SE.mis_0.bam")

# sequence path
plasmid_sequence_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/plasmid_sequences/pCag-EGFP_Lin28IR.fasta"

######################################################## READ DATA
# reads mapped with STAR
star_reads <- #
  GenomicAlignments::readGAlignmentsList(file = star_sample_path, use.names = TRUE, 
                                         param = ScanBamParam(flag = scanBamFlag(isMinusStrand = F), 
                                                              what = "seq")) %>% 
  unlist(.) %>% 
  .[str_detect(cigar(.), "21M|22M|23M")] %>% 
  .[seqnames(.) == "pCag-EGFP_Lin28IR"]

# load reads mapped with shrimp
shrimp_reads <- 
  GenomicAlignments::readGAlignmentsList(file = shrimp_sample_path, use.names = TRUE, 
                                         param = ScanBamParam(flag = scanBamFlag(isMinusStrand = F))) %>% 
  unlist(.) %>% 
  .[str_detect(cigar(.), "21M|22M|23M")] %>% 
  .[seqnames(.) == "pCag-EGFP_Lin28IR"]

# read plasmid sequence 
plasmid_sequence <- Biostrings::readDNAStringSet(plasmid_sequence_path)

######################################################## MAIN CODE
# get read which are missing in STAR
shrimp_only <- shrimp_reads[names(shrimp_reads) %in% setdiff(names(shrimp_reads), names(star_reads))]

# get one read
read_seq <- "TTTGCCAAGTGGCTGGGCTAG"
Biostrings::vmatchPattern(read_seq, plasmid_sequence)
reverseComplement(Biostrings::substr(plasmid_sequence, 3545, 3565))

# shrimp_both <- shrimp_reads[names(shrimp_reads) %in% intersect(names(shrimp_reads), names(star_reads))]
# star_both <- star_reads[names(star_reads) %in% intersect(names(shrimp_reads), names(star_reads))]

#### plot
# reads to plot
bam_gr <- shrimp_only

# get length of chromosome on which is feature located
seq_length <- seqlengths(bam_gr)["pCag-EGFP_Lin28IR"]

# get coverage
coverage_df <- 
  bam_gr %>% 
  coverage(.) %>% 
  .[unique(seqnames(bam_gr))] %>% 
  as(., "IntegerList") %>% 
  unlist(.) %>% 
  unname(.)

if(length(coverage_df) == 0){
  coverage_df <- tibble(pos = 1:seq_length, 
                        coverage = 0)
}else{
  coverage_df <- tibble(pos = 1:seq_length, 
                        coverage = coverage_df)
}

# set position to 0
coverage_df %<>% 
  dplyr::filter((pos >= 1) & (pos <= seq_length)) %>% 
  dplyr::mutate(pos = 1:nrow(.), 
                strand = "plus")

# plot
coverage_plot <-
  ggplot() +
  geom_rect(data = coverage_df, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = coverage, fill = strand)) +
  # coord_cartesian(ylim = c(-200, 200)) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(plot = coverage_plot,
       filename = file.path(outpath, "coverage_plots/whole_plasmids/only_plasmids/compare_aligners", "test.png"),
       width = 15,
       height = 10)


