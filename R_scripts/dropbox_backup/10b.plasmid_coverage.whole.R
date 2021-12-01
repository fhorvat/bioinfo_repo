### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/IR_expression/coverage_plots/whole_plasmid")

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
coverageIR <- function(bam_path, which_gr, on_minus_strand = NA){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = on_minus_strand))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
  bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
  
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

# path to plasmid coordinates
plasmid_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/plasmid_sequences/all_plasmids.bed"

# get path of inverted repeat coordinates from mRNA
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_mRNAs.plasmid_coordinates/IR_mRNA.plasmid_coordinates.clean.csv"

######################################################## READ DATA
# read coordinates of inverted repeats
plasmid_gr <- rtracklayer::import.bed(con = plasmid_path)

# read coordinates of inverted repeats from mRNA
ir_gr <- readr::read_csv(file = ir_path)

######################################################## MAIN CODE
# prepare mRNA IR coordinates on plasmids
ir_gr %<>% GenomicRanges::GRanges(.)
mcols(ir_gr) <- mcols(ir_gr)[, c("arm")]
names(mcols(ir_gr)) <- "gene_id"
mcols(ir_gr)$gene_biotype <- as.character(seqnames(ir_gr))

### get whole plasmid coverage for all samples in each experiment 
# experiment paths
experiment_paths <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" %>% 
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
bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)

# library size path
library_size_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/library_size", str_c("library_size.counts.", experiment, ".csv"))

# read library size table
library_size_df <- 
  readr::read_csv(library_size_path) %>% 
  tidyr::gather(mismatch, library_size, -sample_id) %>% 
  tidyr::separate(mismatch, c("mismatch", "read_group"), sep = "\\.") %>% 
  tidyr::unite(sample_id, sample_id, mismatch, sep = ".") %>% 
  tidyr::spread(read_group, library_size) %>% 
  dplyr::select(sample_id, library_size = `all_reads`)

# loop through sample in experiment 
purrr::map(bam_paths, function(bam_path){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam")
  
  # get plasmid name
  plasmid_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")
  
  # print bam name
  cat(bam_name, "\n")
  
  # subset plasmid ranges
  plasmid_filt <- plasmid_gr[seqnames(plasmid_gr) == plasmid_name]
  
  # subset IR ranges
  ir_filt <-
    ir_gr %>% 
    as.data.frame(.) %>% 
    as.tibble(.) %>% 
    dplyr::filter(seqnames == plasmid_name)
  
  # get coverage of plus and minus strand in one data.frame, normalize for library size
  plot_df <- 
    rbind(coverageIR(bam_path = bam_path, which_gr = plasmid_filt, on_minus_strand = F) %>% 
            dplyr::mutate(strand = "plus"), 
          coverageIR(bam_path = bam_path, which_gr = plasmid_filt, on_minus_strand = T) %>% 
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
    # geom_rect(data = ir_filt, aes(xmin = start, xmax = end, ymin = min(plot_df$rpm), ymax = max(plot_df$rpm)), fill = "red", alpha = 0.1) +
    geom_hline(yintercept = 0, color = "black") +
    # coord_cartesian(xlim = c(min(ir_filt$start), max(ir_filt$end)))+
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save
  ggsave(plot = coverage_plot,
         filename = file.path(outpath, str_c(bam_name, "whole.coverage.RPM.png", sep = ".")),
         width = 15,
         height = 10)
  
  # return
  return(bam_name)
  
})




