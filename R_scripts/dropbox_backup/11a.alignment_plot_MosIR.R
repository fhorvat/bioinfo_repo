### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/phasing/alignment_plots")

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
library(ggbio)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

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
# # ESC bams
# esc_bams <- c("s_ES_Dcr_Rsc_tran_r2.SE.mis_0", "s_ES_Dcr_Rsc_utran_r2.SE.mis_0",
#               "s_ES_DcrSom_Rsc_tran_r1.SE.mis_0", "s_ES_DcrSom_Rsc_utran_r1.SE.mis_0", 
#               "s_ES_WT_tran_r1.SE.mis_0")

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
  
  # loop through sample in experiment 
  for(bam_path in bam_paths){
    
    # get bam name
    bam_name <- basename(bam_path) %>% str_remove(., ".bam")
    
    # print bam name
    cat(bam_name, "\n")
    
    # read bam
    bam_gr <- 
      GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                             use.names = TRUE, 
                                             param = ScanBamParam(which = mosir_gr, 
                                                                  flag = scanBamFlag(isMinusStrand = F))) %>% 
      unlist(.)
    
    # take only reads between 21-23nt
    bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
    
    # add overlap with arm
    values(bam_gr)$arm <- subjectHits(findOverlaps(bam_gr, mosir_gr))
    
    # normalize start of MosIR
    bam_norm <- 
      bam_gr %>% 
      GRanges(.) %>% 
      as.data.frame(.) %>%
      tibble::as.tibble(.) %>%
      dplyr::select(seqnames, start, end, arm) %>% 
      dplyr::mutate(start_norm = ifelse(test = (arm == 1), 
                                        yes = (start - start(mosir_gr[1]) + 1), 
                                        no = (end(mosir_gr[2]) - end + 1)), 
                    end_norm = ifelse(test = (arm == 1), 
                                      yes = (end - start(mosir_gr[1]) + 1),
                                      no = (end(mosir_gr[2]) - start + 1))) %>% 
      dplyr::select(seqnames, start = start_norm, end = end_norm, arm) %>% 
      dplyr::filter(start > 0)
    
    ### collapsed reads
    # get GRanges
    bam_collapsed <- 
      bam_norm %>%
      dplyr::group_by(seqnames, start, end, arm) %>%
      dplyr::summarise(count = n()) %>%
      GRanges(.)
    
    # MosIR arm 1
    reads_plot_arm1 <- 
      ggbio::autoplot(bam_collapsed[bam_collapsed$arm == 1], aes(fill = count, color = count), size = 0.35) +
      scale_fill_continuous(type = "viridis", limits = c(1, max(bam_collapsed$count))) +
      scale_color_viridis_c() +
      guides(color = F)
    
    # MosIR arm 2
    reads_plot_arm2 <- 
      ggbio::autoplot(bam_collapsed[bam_collapsed$arm == 2], aes(fill = count, color = count), size = 0.35) +
      scale_fill_continuous(type = "viridis", limits = c(0, max(bam_collapsed$count))) +
      scale_color_viridis_c() +
      scale_y_reverse(breaks = NULL, minor_breaks = NULL, labels = NULL) +
      guides(fill = F, color = F)
    
    # create and save ggbio plot
    ggbio::tracks(`mosIR arm 1` = reads_plot_arm1,
                  `mosIR arm 2` = reads_plot_arm2,
                  heights = c(1, 1), 
                  theme = theme_clear(), 
                  label.text.cex = 0.75) +
      ggsave(filename = file.path(outpath, str_c(experiment, bam_name, "collapsed", "MosIR_alignments.png", sep = ".")), width = 20, height = 10)
    
    
    # ### individual reads
    # # get GRanges
    # bam_individual <- 
    #   bam_norm %>%
    #   GRanges(.)
    # 
    # # individual reads - MosIR arm 1
    # reads_plot_arm1 <- ggbio::autoplot(bam_individual[bam_individual$arm == 1])
    # 
    # # individual reads - MosIR arm 2
    # reads_plot_arm2 <- 
    #   ggbio::autoplot(bam_individual[bam_individual$arm == 2]) +
    #   scale_y_reverse(breaks = NULL, minor_breaks = NULL, labels = NULL)
    # 
    # # create and save ggbio plot
    # ggbio::tracks(`mosIR arm 1` = reads_plot_arm1,
    #               `mosIR arm 2` = reads_plot_arm2,
    #               heights = c(1, 1), 
    #               theme = theme_clear(), 
    #               label.text.cex = 0.75) +
    #   ggsave(filename = file.path(outpath, str_c(experiment, bam_name, "individual", "MosIR_alignments.png", sep = ".")), width = 20, height = 10)
    
  }
  
  cat("\n")
  
}



