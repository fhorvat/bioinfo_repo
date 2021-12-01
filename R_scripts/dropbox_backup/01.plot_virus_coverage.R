### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2020_Apr.brain_spleen/Analysis/virus_genome_coverage")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get coverage from bam file
coverageSequence <- function(bam_path, which_gr, minus.strand = F){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = minus.strand))) %>% 
    unlist(.)
  
  # # take only reads between 21-23nt
  # bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
  
  # get length of chromosome on which is feature located
  seq_length <- seqlengths(bam_gr)[as.character(seqnames(which_gr))]
  
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
    dplyr::filter((pos >= start(which_gr)) & (pos <= end(which_gr))) %>% 
    dplyr::mutate(pos = 1:nrow(.)) %>% 
    dplyr::mutate(strand = ifelse(minus.strand, "-", "+"))
  
  # return 
  return(coverage_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2020_Apr.brain_spleen/Data/Mapped/STAR_mm10"

######################################################## READ DATA

######################################################## MAIN CODE
# set whole plasmid coordinates
viral_gr <- GenomicRanges::GRanges(seqnames = "NC_001672.1",
                                   ranges = IRanges(start = 1, end = 11141),
                                   gene_id = "NC_001672.1")

### get viral genome counts for all samples, normalize, plot
# list bams 
bam_paths <- list.files(path = mapped_path, pattern = "s_brain.*\\.bam$", full.names = T)

# library size path
library_size_path <- file.path(mapped_path, "library_sizes.txt")

# read library size df
library_size_df <- 
  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) %>% 
  tibble::as_tibble(.)

# loop through sample in experiment 
for(bam_path in bam_paths){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., "\\.bam")
  
  # print bam name
  cat(bam_name, "\n")
  
  # get coverage of both mosIR arms in one data.frame, normalize for library size
  plot_df <-
    dplyr::bind_rows(coverageSequence(bam_path = bam_path, which_gr = viral_gr, minus.strand = T), 
                     coverageSequence(bam_path = bam_path, which_gr = viral_gr, minus.strand = F)) %>% 
    dplyr::mutate(sample_id = bam_name) %>%
    dplyr::left_join(., library_size_df, by = "sample_id") %>%
    dplyr::mutate(library_size = round((library_size / 1e6), 4),
                  rpm = coverage / library_size) %>% 
    dplyr::mutate(rpm = ifelse(strand == "+", rpm, -rpm))
  
  # get plot limits
  plot_lim <- plot_df$rpm %>% abs %>% max %>% ceiling
  
  # plot
  coverage_plot <-
    ggplot() +
    geom_rect(data = plot_df, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
    geom_hline(yintercept = 0, color = "black") +
    # coord_cartesian(ylim = c(-plot_lim, plot_lim)) +
    guides(fill = FALSE) +
    ggtitle(bam_name) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save
  ggsave(plot = coverage_plot,
         filename = file.path(outpath, str_c(bam_name, "viral_genome_coverage.RPM.png", sep = ".")),
         width = 15,
         height = 10)
  
}
