### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/Mos_coverage")

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
library(Biostrings)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# get coverage from bam file (with 0, 1 or 2 mismatches)
bamCoverage <- function(bam_path, isMinusStrand = NA){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".SE.genome.Aligned.sortedByCoord.out.bam")
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(bam_path, use.names = TRUE, param = ScanBamParam(tag = c("nM", "NH"), flag = scanBamFlag(isMinusStrand = isMinusStrand))) %>% 
    unlist(.)
  
  # filter bam by number of mismatches
  coverage_mismatch <- purrr::map(c(1, 2, 3), function(mismatch){
    
    # filter
    bam_gr <- bam_gr[mcols(bam_gr)$nM < mismatch]
    
    # get coverage
    coverage_bam <- 
      bam_gr %>% 
      coverage(.) %>% 
      as(., "IntegerList") %>% 
      unlist(.) %>% 
      unname(.) %>% 
      tibble(.)
    
    # return 
    return(coverage_bam)
    
  }) %>% 
    dplyr::bind_cols(.) %>% 
    magrittr::set_colnames(c("mismatch_0", "mismatch_1", "mismatch_2"))
  
  
  # add to tibble
  coverage_df <- 
    coverage_mismatch %>% 
    dplyr::mutate(pos = 1:nrow(.), 
                  sample_id = bam_name) %>% 
    dplyr::select(pos, everything())
  
  # return 
  return(coverage_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get bam file path, name, experiment
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_pCMV_MosIR_EGFP/filtered_21to23nt"

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*SE.21to23nt.bam$", full.names = T)

# library size path
library_size_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10_new/filtered_21to23nt/library_sizes.21to23nt.txt"

# documentation path
documentation_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Documentation")

# sample table path
sample_table_path <- list.files(path = documentation_path, ".*sampleTable.csv", full.names = T)

# pCMV MosIR EGFP sequence path
vector_seq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Documentation/pCMV_MosIR_EGFP/pCMV_MosIR_EGFP.fasta"
  
######################################################## READ DATA
# read library size table
library_size_df <- 
  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, ".SE.21to23nt"), 
                library_size = round((library_size / 1e6), 3))

# read sample table
sample_table <- 
  readr::read_csv(file = sample_table_path) %>% 
  dplyr::select(sample_id, genotype, transfection) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, ".SE"))

# read vector sequence
vector_seq <- 
  readDNAStringSet(vector_seq_path) %>% 
  unlist(.)

######################################################## MAIN CODE
# get sequence of first and second arm
vector_seq_1 <- 
  Biostrings::substring(vector_seq, 1800, 2319) %>% 
  as.character(.)
vector_seq_2 <- 
  Biostrings::substring(vector_seq, 2418, 2937) %>% 
  reverseComplement(.) %>% 
  as.character(.)

# get coverage of all bam files in one data.frame
coverage_counts <- 
  purrr::map(bam_paths, bamCoverage, isMinusStrand = F) %>% 
  dplyr::bind_rows(.)

# normalize coverage to CPM
coverage_cpm <- 
  coverage_counts %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, ".SE.21to23nt.bam")) %>% 
  dplyr::left_join(., library_size_df, by = "sample_id") %>% 
  dplyr::mutate_at(.vars = vars(starts_with("mismatch")), .funs = funs(round((. / library_size), 3)))

# get max CPM per genotype/transfection combination
scale_CPM <- 
  coverage_cpm %>% 
  dplyr::left_join(., sample_table, by = "sample_id") %>% 
  dplyr::group_by(genotype, transfection) %>% 
  dplyr::summarise_at(.vars = vars(starts_with("mismatch")), .funs = funs(max(.))) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::right_join(sample_table)

### plot coverage 
# plot faceted  
for(mismatch in c(0, 1, 2)){
  
  # data.frame for plot
  plot_df <- 
    coverage_cpm %>% 
    dplyr::select(x = pos, y = str_c("mismatch_", mismatch), sample_id)
  
  # plot 
  coverage_plot <- 
    ggplot(plot_df, aes(x, y, sample_id)) +
    geom_rect(aes(xmin = x, xmax = x + 1, ymin = 0, ymax = y), fill = "black", color = "black") +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # plot with free_y axis 
  coverage_plot +
    facet_grid(sample_id ~ ., scales = "free_y") + 
    ggsave(filename = file.path(outpath, str_c("facet.pCMV_MosIR_EGFP.CPM.coverage.", mismatch, "_mismatches_allowed.pdf")), width = 10, height = 30)
  
  # # plot with scaled y-axis 
  # coverage_plot +
  #   facet_grid(sample_id ~ .) + 
  #   ggsave(filename = file.path(outpath, str_c("pCMV_MosIR_EGFP.coverage.", mismatch, "_mismatches_allowed.pdf")), width = 10, height = 30)
  
}

# plot separate samples
for(sample in unique(coverage_cpm$sample_id)){
  
  for(mismatch in c(0, 1, 2)){
    
    # data.frame for plot
    plot_df <- 
      coverage_cpm %>% 
      dplyr::filter(sample_id == sample) %>% 
      dplyr::select(x = pos, y = str_c("mismatch_", mismatch), sample_id)
    
    # max Y scale value
    max_y <- 
      scale_CPM %>% 
      dplyr::filter(sample_id == sample) %>% 
      # dplyr::select(max_y = str_c("mismatch_", mismatch)) %$% 
      dplyr::select(max_y = str_c("mismatch_", 2)) %$% 
      max_y
    
    # plot 
    coverage_plot <- 
      ggplot(plot_df, aes(x, y, sample_id)) +
      geom_rect(aes(xmin = x, xmax = x + 1, ymin = 0, ymax = y), fill = "black", color = "black") +
      scale_y_continuous(limits = c(0, plyr::round_any(max_y, 10, f = ceiling)), 
                         breaks = seq(0, plyr::round_any(max_y, 10, f = ceiling), 10)) +
      theme_bw() +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            axis.text.x = element_text(size = 10), 
            axis.text.y = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ggsave(filename = file.path(outpath, str_c(sample, ".pCMV_MosIR_EGFP", ".CPM.coverage.", mismatch, "_mismatches_allowed.pdf")), width = 10, height = 2)
    
  }
  
}

# plot faceted with limited scale
for(mismatch in c(0, 1)){
  
  # data.frame for plot
  plot_df <- 
    coverage_cpm %>% 
    dplyr::select(x = pos, y = str_c("mismatch_", mismatch), sample_id)
  
  # plot 
  coverage_plot <- 
    ggplot(plot_df, aes(x, y, sample_id)) +
    geom_rect(aes(xmin = x, xmax = x + 1, ymin = 0, ymax = y), fill = "black", color = "black") +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_grid(sample_id ~ .)
  
  # plot with linear limited y-axis 
  coverage_plot +
    coord_cartesian(ylim = c(0, 2000)) +
    ggsave(filename = file.path(outpath, str_c("facet.pCMV_MosIR_EGFP.CPM.coverage.sense.", mismatch, "_mismatches_allowed.linear_2000.pdf")), width = 10, height = 30)
  
  # plot with logaritmic limited y-axis 
  coverage_plot +
    scale_y_log10() + 
    coord_cartesian(ylim = c(0.1, 2000)) +
    ggsave(filename = file.path(outpath, str_c("facet.pCMV_MosIR_EGFP.CPM.coverage.sense.", mismatch, "_mismatches_allowed.log10_2000.pdf")), width = 10, height = 30)
  
  
}
