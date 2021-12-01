### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression/limited_plots")

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
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# get coverage from bam file
coverageMosIR <- function(bam_path, which_gr){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = F))) %>% 
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

# MosIR .bed path
mosir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

######################################################## READ DATA
# read bed with coordinates
mosir_gr <- rtracklayer::import.bed(con = mosir_path)
mosir_gr <- mosir_gr[mosir_gr$name != "EGFP"]

######################################################## MAIN CODE
# # set whole plasmid coordinates
# plasmid_gr <- GenomicRanges::GRanges(seqnames = "pCAG-EGFP_MosIR", 
#                                      ranges = IRanges(start = 1, end = 10000), 
#                                      gene_id = "pCAG-EGFP_MosIR")

### get MosIR counts for all samples in each experiment 
# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/T3T_mESC_MosIR.2016/Data/Mapped/STAR_mm10",
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Mapped/STAR_mm10.pCAG_EGFP_MosIR") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  # print experiment
  cat("Plotting experiment:", experiment, "\n")
  
  # get mapped path 
  mapped_path <- experiment_paths[experiment]
  
  # list bams 
  bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)
  
  # library size path
  library_size_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/MosIR_expression", 
                                 str_c("mosIR.", experiment, ".counts_summary.xlsx"))
  
  # read library size df
  library_size_df <- 
    xlsx::read.xlsx(file = library_size_path, sheetName = "library_sizes") %>% 
    tibble::as.tibble(.) %>% 
    tidyr::gather(library_type, library_size, -sample_id) %>% 
    tidyr::separate(library_type, into = c("mis", "library_type"), sep = "\\.") %>% 
    tidyr::unite(sample_id, sample_id, mis, sep = ".") %>% 
    tidyr::spread(library_type, library_size)
  
  # loop through sample in experiment 
  for(bam_path in bam_paths){
    
    # get bam name
    bam_name <- basename(bam_path) %>% str_remove(., ".bam")
    
    # print bam name
    cat(bam_name, "\n")
    
    # get coverage of both mosIR arms in one data.frame, normalize for library size
    plot_df <-
      purrr::map(.x = 1:2, function(x){
        
        # get coverage
        coverage_df <-
          coverageMosIR(bam_path = bam_path, which_gr = mosir_gr[x]) %>%
          dplyr::mutate(arm = mcols(mosir_gr[x])$name)
        
        # reverse coverage for second MosIR arm
        if(x == 2){
          
          coverage_df %<>%
            dplyr::mutate(pos = rev(pos) + 1,
                          coverage = - coverage)
          
        }
        
        return(coverage_df)
        
      }) %>%
      dplyr::bind_rows(.) %>%
      dplyr::mutate(sample_id = bam_name) %>%
      dplyr::left_join(., library_size_df %>% dplyr::select(sample_id, library_size = `21to23nt_reads`), by = "sample_id") %>%
      dplyr::mutate(library_size = round((library_size / 1e6), 4),
                    rpm = coverage / library_size)
    
    # plot
    coverage_plot <-
      ggplot() +
      geom_rect(data = plot_df, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = arm)) +
      geom_hline(yintercept = 0, color = "black") +
      coord_cartesian(ylim = c(-1500, 1500)) + 
      guides(fill = FALSE) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # save
    ggsave(plot = coverage_plot,
           filename = file.path(outpath, experiment, str_c(experiment, bam_name, "MosIR_coverage.RPM.lim.png", sep = ".")),
           width = 15,
           height = 10)
    
  }
  
  cat("\n")
  
}

