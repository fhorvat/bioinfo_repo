### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/phasing")

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
library(ggiraphExtra)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS
# get coverage from bam file
phasingRegister <- function(bam_path, which_gr){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = F))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
  bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
  
  # get starting positions of reads, normalize to start of mosIR arm, get mod 21
  if(values(which_gr)$name == "MosIR_arm1"){
    
    # if read is on MosIR arm 1 just set start of arm 1 to 1
    values(bam_gr)$start_norm <- start(bam_gr) - start(which_gr) + 1
    values(bam_gr)$end_norm <- end(bam_gr) - start(which_gr) + 1
    
  }else{
    if(values(which_gr)$name == "MosIR_arm2"){
      
      # if read is on MosIR arm 2 then start is really an end of mosIR arm 2
      values(bam_gr)$start_norm <- end(which_gr) - end(bam_gr) + 1
      values(bam_gr)$end_norm <- end(which_gr) - start(bam_gr) + 1
      
    }else{
      
      # if non true then input is wrong
      stop("Something went wrong!")
      
    }
    
  }
  
  # get mod 21 of normalized read starts
  read_start_mod <- (values(bam_gr)$start_norm %% 22) + 1
  
  # table for plot
  read_register_df <- 
    tibble(read_register = read_start_mod) %>% 
    dplyr::group_by(read_register) %>% 
    dplyr::summarize(count = n()) %>% 
    dplyr::mutate(percent = (count / sum(count)) * 100) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-count) %>% 
    dplyr::right_join(tibble(read_register = 1:22), by = "read_register") %>% 
    dplyr::mutate(percent = replace(percent, is.na(percent), 0)) %>% 
    tidyr::spread(read_register, percent)
  
  # return 
  return(read_register_df)
  
}

roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
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
start(mosir_gr) <- start(mosir_gr) - 1

######################################################## MAIN CODE
# ESC bams
esc_bams <- c("s_ES_Dcr_Rsc_tran_r2.SE.mis_0", "s_ES_Dcr_Rsc_utran_r2.SE.mis_0",
              "s_ES_DcrSom_Rsc_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_utran_r1.SE.mis_0", 
              "s_ES_WT_tran_r1.SE.mis_0")

# experiment paths
experiment_paths <- 
  c("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/T3T_DcrTrans_2011/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/SOLiD/ES_DcrTrans_2012/Data/Mapped/Shrimp_mm10", 
    "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/T3T_mESC_MosIR.2016/Data/Mapped/STAR_mm10") %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_paths)

# loop through experiments
for(experiment in experiment_list){
  
  experiment <- "T3T_mESC_MosIR.2016"
    
  # print experiment
  cat("Plotting experiment:", experiment, "\n")
  
  # get mapped path 
  mapped_path <- experiment_paths[experiment]
  
  # list bams 
  bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)
  
  # filter bams in ESC
  # if(experiment == "ES_DcrTrans_2012"){
  #   bam_paths <- bam_paths[str_detect(bam_paths, str_c(esc_bams, collapse = "|"))]
  # }
  
  # loop through sample in experiment 
  plot_facet_df <- purrr::map(.x = bam_paths, function(bam_path){
    
    # get bam name
    bam_name <- basename(bam_path) %>% str_remove(., ".bam")
    
    # print bam name
    cat(bam_name, "\n")
    
    # get coverage of both mosIR arms in one data.frame, normalize for library size
    plot_df <-
      purrr::map(.x = 1:2, function(x){
        
        # get coverage
        phasing_df <-
          phasingRegister(bam_path = bam_path, which_gr = mosir_gr[x]) %>%
          dplyr::mutate(arm = mcols(mosir_gr[x])$name)
        
        return(phasing_df)
        
      }) %>%
      dplyr::bind_rows(.) %>% 
      dplyr::mutate(sample_id = bam_name) %>% 
      dplyr::select(sample_id, arm, everything())
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # create and save radar plot
  ggRadar(data = plot_facet_df, aes(color = arm, facet = sample_id), 
          rescale = F, interactive = F, size = 0.1, 
          ylim = c(0, 30)) +
    facet_wrap(~ sample_id) +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10), 
          panel.grid.major = element_line(size = 1, color = "gray90")) +
    ggsave(filename = file.path(outpath, str_c(experiment, "MosIR_phase.facet.png", sep = ".")), width = 20, height = 20)
  
  # new line
  cat("\n")
  
  
  ### plot as separate samples
  # split data.frame to list
  plot_list <- 
    plot_facet_df %>% 
    split(., plot_facet_df$sample_id)
  
  # plot in a loop
  invisible(purrr::map(names(plot_list), function(bam_name){
    
    # print bam name
    cat(bam_name, "\n")
    
    # get plot data.frame
    plot_df <- plot_list[[bam_name]]
      
    # create and save radar plot
    ggRadar(data = plot_df, aes(color = arm), 
            rescale = F, interactive = F, size = 0.1, 
            ylim = c(0, 30)) +
      theme_bw() +
      theme(legend.position = "bottom", 
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10), 
            panel.grid.major = element_line(size = 1, color = "gray90")) +
      ggsave(filename = file.path(outpath, str_c(experiment, bam_name, "MosIR_phase.png", sep = ".")), width = 10, height = 10)
    
    # return
    return(bam_name)
    
  }))
  
}



