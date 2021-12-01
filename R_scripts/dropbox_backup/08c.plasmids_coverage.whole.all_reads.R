### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/IR_expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(xlsx)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
### scans chromosomes in bam file
chrFinder <- function(bam.path, filter = FALSE, output = "data.frame"){
  
  # scan bam
  s <- scanBamHeader(bam.path)
  st <- s[[1]]$text
  st <- do.call(rbind, st[names(st) == "@SQ"])
  st[, 1] <- str_replace(st[,1], "SN:", "")
  st[, 2] <- str_replace(st[,2], "LN:", "")
  
  # filter
  if(filter == TRUE){
    st <- st[!str_detect(st[, 1], "random")]
  }
  
  # output 
  if(output == 'data.frame'){
    
    vst <- data.frame(chr = st[, 1], chrlen = as.numeric(st[, 2]), stringsAsFactors = F)
    
  }else{
    
    vst <- as.numeric(st[, 2])
    names(vst) <- st[, 1]
    
  }
  
  # return
  return(vst)
  
}

# get coverage from bam file
coverageIR <- function(bam_path, which_gr, on_minus_strand = NA){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = on_minus_strand))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
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

# get .BED with inverted repeat coordinates
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/plasmid_sequences/all_plasmids.bed"

# sample bam path
sample_bam_path <- list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.all_plasmids", 
                              pattern = "*.bam$", full.names = T)[1]

######################################################## READ DATA
# read coordinates of inverted repeats
ir_gr <- rtracklayer::import.bed(con = ir_path)

######################################################## MAIN CODE
# get coordinates of whole plasmids
plasmid_gr <- 
  chrFinder(sample_bam_path) %>% 
  dplyr::filter(str_detect(chr, "^pCag-EGFP|^pU6")) %>% 
  dplyr::mutate(start = 1, 
                gene_id = chr) %>% 
  dplyr::select(seqnames = chr, start, end = chrlen, gene_id) %>% 
  GenomicRanges::GRanges(.)

### get MosIR counts for all samples in each experiment 
# experiment paths
experiment_paths <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.all_plasmids" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
experiment <- experiment_list[1]

# print experiment
cat("Plotting experiment:", experiment, "\n")

# get mapped path 
mapped_path <- experiment_paths[experiment]

# list bams 
bam_paths <- list.files(path = file.path(mapped_path, "3_original_mapping"), pattern = ".*bam$", full.names = T)

# library size path
library_size_path <- file.path(mapped_path, "4_library_size/library_sizes.txt")

# read library size df
library_size_df <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

# loop through sample in experiment 
for(bam_path in bam_paths){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam")
  
  # print bam name
  cat(str_c(bam_name, ": "))
  
  for(n in 1:length(plasmid_gr)){
    
    # subset plasmid ranges
    plasmid <- plasmid_gr[n]
    plasmid_name <- seqnames(plasmid)
    cat(as.character(plasmid_name), "")
    
    # get coverage of plus and minus strand in one data.frame, normalize for library size
    plot_df <- 
      rbind(coverageIR(bam_path = bam_path, which_gr = plasmid, on_minus_strand = F) %>% 
              dplyr::mutate(strand = "plus"), 
            coverageIR(bam_path = bam_path, which_gr = plasmid, on_minus_strand = T) %>% 
              dplyr::mutate(strand = "minus", 
                            coverage = - coverage)) %>% 
      dplyr::mutate(sample_id = bam_name) %>%
      dplyr::left_join(., library_size_df, by = "sample_id") %>%
      dplyr::mutate(library_size = round((library_size / 1e6), 4),
                    rpm = coverage / library_size)
    
    # plot
    coverage_plot <-
      ggplot() +
      geom_rect(data = plot_df, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
      geom_hline(yintercept = 0, color = "black") +
      # coord_cartesian(ylim = c(-1500, 1500)) +
      guides(fill = FALSE) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # save
    ggsave(plot = coverage_plot,
           filename = file.path(outpath, "coverage_plots/whole_plasmids/all_reads", str_c(bam_name, plasmid_name, "whole.coverage.RPM.png", sep = ".")),
           width = 15,
           height = 10)
    
  }
  
  cat("\n")
  
}


