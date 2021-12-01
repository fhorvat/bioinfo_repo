### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb/Data/Mapped/STAR_DicerLocus")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get coverage from bam file
coverageScaffold <- function(bam_path, which_gr, on_minus_strand){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = on_minus_strand))) %>% 
    unlist(.)
  
  # remove reads with deletions and insertions
  bam_gr <- bam_gr[!str_detect(cigar(bam_gr), "D|I")]
  # bam_gr <- bam_gr[!str_detect(cigar(bam_gr), "S")]
  
  # # get reads longer than 30 nt
  bam_gr <- bam_gr[width(bam_gr) > 20]
  
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
    dplyr::mutate(pos = 1:nrow(.))
  
  # return 
  return(coverage_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# .bed path
bed_path <- file.path(inpath, "STAR_index", "Dicer1.locus.bed")

# mapped path
mapped_path <- inpath

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*\\.bam$", full.names = T)

######################################################## READ DATA
# read bed with coordinates
bed_gr <- rtracklayer::import.bed(con = bed_path)

######################################################## MAIN CODE
### prepare file
bed_tb <- as_tibble(bed_gr) %>% dplyr::rename(locus = seqnames)
bed_gr <- reduce(bed_gr)


### create table with coverage
# loop through samples in experiment 
coverage_tb_full <- purrr::map(bam_paths, function(path){
  
  # get bam name
  bam_name <- basename(path) %>% str_remove(., ".bam")
  
  # print bam name
  cat(bam_name, "\n")
  
  # get coverage of both mosIR arms in one data.frame, normalize for library size
  coverage_tb <- 
    rbind(coverageScaffold(bam_path = path, which_gr = bed_gr[1], on_minus_strand = F) %>% 
            dplyr::mutate(locus = as.character(seqnames(bed_gr[1]))),
          coverageScaffold(bam_path = path, which_gr = bed_gr[2], on_minus_strand = F) %>% 
            dplyr::mutate(locus = as.character(seqnames(bed_gr[2])))) %>% 
    dplyr::mutate(sample_id = bam_name)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# create table for plot
plot_tb <- 
  coverage_tb_full %>% 
  dplyr::mutate(genotype_dicer = str_extract(sample_id, "DicerX_HET|WT"), 
                genotype_transfection = str_extract(sample_id, "TBEV")) %>% 
  dplyr::mutate(genotype_transfection = replace(genotype_transfection, is.na(genotype_transfection), "no_transfection")) %>% 
  dplyr::mutate(genotype_dicer = factor(genotype_dicer, levels = c("DicerX_HET", "WT")), 
                genotype_transfection = factor(genotype_transfection, levels = c("TBEV", "no_transfection"))) %>% 
  dplyr::arrange(genotype_dicer, genotype_transfection) %>% 
  dplyr::mutate(sample_name = str_c(str_remove(genotype_dicer, "Dicer "), str_extract(sample_id, "[5-9]M"), sep = " "))

# get y-limit
ylimit <- plot_tb$coverage %>% abs(.) %>% max(.) %>% ceiling(.)

### plot whole virus
# plot
coverage_plot <-
  ggplot() +
  geom_rect(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = coverage), fill = "black") +
  geom_rect(data = bed_tb, aes(xmin = start, xmax = end, ymin = -1, ymax = 0, fill = name)) +
  scale_fill_manual(values = c("exon2" = "#009900", "HAtag" = "#ff0000", "exon7" = "#0000ff", "exon3" = "orange")) +
  facet_grid(sample_name ~ locus, scales = "free_y") + 
  # guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 10), 
        strip.text.y = element_text(size = 15))

# save
ggsave(plot = coverage_plot,
       filename = file.path(outpath, str_c("Dicer_locus_coverage.20plus_reads.png", sep = ".")),
       width = 20,
       height = 15)
 