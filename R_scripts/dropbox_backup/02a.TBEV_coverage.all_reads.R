### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb/Analysis/TBEV_coverage")

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
    dplyr::mutate(pos = 1:nrow(.))
  
  # return 
  return(coverage_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# MosIR .bed path
bed_path <- file.path(inpath, "NC_001672.1.TBEV.bed")

# experiment paths
experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.TBEV.small_RNAseq.2021_Feb"

# mapped path
mapped_path <- file.path(experiment_path, "Data/Mapped/STAR_mm10")

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*\\.bam$", full.names = T)

# library size path
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

######################################################## READ DATA
# read bed with coordinates
bed_gr <- rtracklayer::import.bed(con = bed_path)

# read library size df
library_size_df <-  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

######################################################## MAIN CODE
### prepare files
# clean library size
library_size_tidy <- 
  library_size_df %>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt$")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"))


### create table with coverage
# loop through samples in experiment 
coverage_tb_full <- purrr::map(bam_paths, function(path){
  
  # get bam name
  bam_name <- basename(path) %>% str_remove(., ".bam")
  
  # print bam name
  cat(bam_name, "\n")
  
  # get coverage of both mosIR arms in one data.frame, normalize for library size
  coverage_tb <- 
    rbind(coverageScaffold(bam_path = path, which_gr = bed_gr, on_minus_strand = F) %>% 
            dplyr::mutate(strand = "plus"), 
          coverageScaffold(bam_path = path, which_gr = bed_gr, on_minus_strand = T) %>% 
            dplyr::mutate(strand = "minus", 
                          coverage = - coverage)) %>% 
    dplyr::mutate(sample_id = bam_name) %>%
    dplyr::left_join(., library_size_df, by = "sample_id") %>%
    dplyr::mutate(library_size = round((library_size / 1e6), 4),
                  rpm = coverage / library_size)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# create table for plot
plot_tb <- 
  coverage_tb_full %>% 
  dplyr::filter(str_detect(sample_id, "\\.SE$")) %>% 
  dplyr::mutate(genotype_dicer = str_extract(sample_id, "DicerX_HET|WT"), 
                genotype_transefection = str_extract(sample_id, "TBEV")) %>% 
  dplyr::mutate(genotype_transefection = replace(genotype_transefection, is.na(genotype_transefection), "no_transfection")) %>% 
  dplyr::mutate(genotype_dicer = factor(genotype_dicer, levels = c("DicerX_HET", "WT")), 
                genotype_transefection = factor(genotype_transefection, levels = c("TBEV", "no_transfection"))) %>% 
  dplyr::arrange(genotype_dicer, genotype_transefection)

# get y-limit
ylimit <- plot_tb$rpm %>% abs(.) %>% max(.) %>% ceiling(.)

### plot whole virus
# plot
coverage_plot <-
  ggplot() +
  geom_rect(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(vars(sample_id), nrow = length(unique(plot_tb$sample_id))) + 
  scale_y_continuous(limits = c(-ylimit, ylimit)) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(plot = coverage_plot,
       filename = file.path(outpath, str_c("TBEV_coverage.all_reads.whole_virus.RPM.png", sep = ".")),
       width = 8,
       height = 15)


### plot whole virus - zoom
# plot
coverage_plot <-
  ggplot() +
  geom_rect(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(vars(sample_id), nrow = length(unique(plot_tb$sample_id))) + 
  coord_cartesian(ylim  = c(-0.1, 0.1)) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(plot = coverage_plot,
       filename = file.path(outpath, str_c("TBEV_coverage.all_reads.zoom.whole_virus.RPM.png", sep = ".")),
       width = 8,
       height = 15)


