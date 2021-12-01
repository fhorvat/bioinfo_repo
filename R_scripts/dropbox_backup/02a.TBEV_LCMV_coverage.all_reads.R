### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.mESC.TBEV_LCMV.small_RNAseq.2021_May/Analysis/coverage")

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
  
  # get coverage
  coverage_all <- 
    bam_gr %>% 
    coverage(.) %>% 
    .[unique(seqnames(bam_gr))]
  
  # loop through different seqnames
  coverage_tb <- purrr::map(as.character(seqnames(which_gr)), function(seqname){
    
    # get one sequence coverage
    coverage_df <- 
      coverage_all[[seqname]] %>% 
      as(., "IntegerList") %>% 
      unlist(.) %>% 
      unname(.)
    
    # get length of chromosome on which is feature located
    seq_length <- seqlengths(bam_gr)[seqname]
    
    # create table
    if(length(coverage_df) == 0){
      coverage_df <- tibble(pos = 1:seq_length, 
                            coverage = 0)
    }else{
      coverage_df <- tibble(pos = 1:seq_length, 
                            coverage = coverage_df)
    }
    
    # set position to 0
    coverage_df %<>% 
      dplyr::filter((pos >= start(which_gr[seqnames(which_gr) == seqname])) & (pos <= end(which_gr[seqnames(which_gr) == seqname]))) %>% 
      dplyr::mutate(pos = 1:nrow(.)) %>% 
      dplyr::mutate(seq_names = seqname)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  
  # return 
  return(coverage_tb)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# MosIR .bed path
bed_path <- file.path(inpath, "viral_sequences.bed")

# experiment paths
experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/antiviral_RNAi.Marcos/datasets/mouse_brain.mESC.TBEV_LCMV.small_RNAseq.2021_May"

# mapped path
mapped_path <- file.path(experiment_path, "Data/Mapped/STAR_mm10_TBEV_LCMV")

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*\\SE.bam$", full.names = T)

# library size path
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# sample table path
sample_table_path <- file.path(experiment_path, "Data/Documentation")
sample_table_path <- list.files(sample_table_path, ".*\\.sampleTable.csv", full.names = T)

######################################################## READ DATA
# read bed with coordinates
bed_gr <- rtracklayer::import.bed(con = bed_path)

# read library size df
library_size_df <-  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
### prepare files
# clean library size
library_size_tidy <- 
  library_size_df %>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt$")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"))

# clean sample table
sample_table_tidy <- 
  sample_table %>% 
  dplyr::select(sample_id, tissue, genotype, infection) %>% 
  dplyr::mutate(infection = replace(infection, infection == "LCMV", "LMCV"))
  
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
plot_tb_list <- 
  coverage_tb_full %>% 
  dplyr::left_join(., sample_table_tidy, by = "sample_id") %>% 
  dplyr::mutate(tissue = factor(tissue, levels = c("brain", "RS10_mESC", "RS7_mESC")), 
                genotype = factor(genotype, levels = c("DicerX", "DicerOO_PKRdel", "WT")), 
                infection = factor(infection, levels = c("TBEV", "LMCV"))) %>% 
  dplyr::arrange(tissue, genotype, infection)

# loop through different seqnames
purrr::map(unique(coverage_tb_full$seq_names), function(seq_name){
  
  # get infection
  true_infection <- str_extract(seq_name, "LMCV|TBEV")
  
  # filter
  plot_tb <- 
    plot_tb_list %>% 
    dplyr::filter(seq_names == seq_name) %>% 
    dplyr::filter(infection == true_infection)
  
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
         filename = file.path(outpath, str_c(seq_name, "coverage", "all_reads.RPM.png", sep = ".")),
         width = 8,
         height = 15)
  
  
  # ### plot whole virus - zoom
  # # plot
  # coverage_plot <-
  #   ggplot() +
  #   geom_rect(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
  #   geom_hline(yintercept = 0, color = "black") +
  #   facet_wrap(vars(sample_id), nrow = length(unique(plot_tb$sample_id))) + 
  #   coord_cartesian(ylim  = c(-0.1, 0.1)) +
  #   guides(fill = FALSE) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(size = 10),
  #         axis.text.y = element_text(size = 10),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # # save
  # ggsave(plot = coverage_plot,
  #        filename = file.path(outpath, str_c("TBEV_coverage.all_reads.zoom.whole_virus.RPM.png", sep = ".")),
  #        width = 8,
  #        height = 15)
  
  # return
  return(seq_name)
  
})



