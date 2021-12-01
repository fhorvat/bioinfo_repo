### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

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
  if(str_detect(values(which_gr)$gene_id, "\\.antisense_arm")){
    
    # if read is on MosIR arm 1 just set start of arm 1 to 1
    values(bam_gr)$start_norm <- start(bam_gr) - start(which_gr) + 1
    values(bam_gr)$end_norm <- end(bam_gr) - start(which_gr) + 1
    
  }else{
    if(str_detect(values(which_gr)$gene_id, "\\.sense_arm")){
      
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

roundUpNice <- function(x, nice = c(1, 2, 3, 4, 5, 6, 8, 10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get path of inverted repeat coordinates from mRNA
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_mRNAs.plasmid_coordinates/IR_mRNA.plasmid_coordinates.clean.csv"

######################################################## READ DATA
# read coordinates of inverted repeats from mRNA
ir_gr <- readr::read_csv(file = ir_path)

######################################################## MAIN CODE
# prepare mRNA IR coordinates on plasmids
ir_gr %<>% GenomicRanges::GRanges(.)
mcols(ir_gr) <- mcols(ir_gr)[, c("arm")]
names(mcols(ir_gr)) <- "gene_id"
mcols(ir_gr)$gene_biotype <- as.character(seqnames(ir_gr))

# experiment paths
experiment_path <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_path)

# loop through experiments
experiment <- experiment_list[1]

# print experiment
cat("Plotting experiment:", experiment, "\n")

# get mapped path 
mapped_path <- experiment_path[experiment]

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)

# loop through sample in experiment 
plot_facet_df <- purrr::map(.x = bam_paths, function(bam_path){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam")
  
  # get IR name
  ir_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")
  
  # subset IR GenomicRanges
  ir_filt <- ir_gr[seqnames(ir_gr) == ir_name]
  
  # print bam name
  cat(bam_name, "\n")
  
  # get coverage of both mosIR arms in one data.frame, normalize for library size
  plot_df <-
    purrr::map(.x = 1:2, function(x){
      
      # get coverage
      phasing_df <-
        phasingRegister(bam_path = bam_path, which_gr = ir_filt[x]) %>%
        dplyr::mutate(arm = mcols(ir_filt[x])$gene_id)
      
      return(phasing_df)
      
    }) %>%
    dplyr::bind_rows(.) %>% 
    dplyr::mutate(sample_id = bam_name) %>% 
    dplyr::select(sample_id, arm, everything())
  
}) %>% 
  dplyr::bind_rows(.)

# # create and save radar plot
# ggRadar(data = plot_facet_df, aes(color = arm, facet = sample_id), 
#         rescale = F, interactive = F, size = 0.1) + 
#   facet_wrap(~ sample_id) +
#   theme_bw() +
#   theme(legend.position = "bottom", 
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10), 
#         panel.grid.major = element_line(size = 1, color = "gray90")) +
#   ggsave(filename = file.path(outpath, experiment, str_c("IR_phase.facet.png", sep = ".")), width = 20, height = 20)


### plot as separate samples
# split data.frame to list
plot_list <- 
  plot_facet_df %>% 
  split(., plot_facet_df$sample_id)

# plot in a loop
invisible(purrr::map(names(plot_list), function(bam_name){
  
  bam_name <- names(plot_list)[1]
  
  # print bam name
  cat(bam_name, "\n")
  
  # get plot data.frame
  plot_df <- plot_list[[bam_name]]
  
  # create and save radar plot
  ggRadar(data = plot_df, mapping = aes(color = arm), 
          rescale = F, interactive = F, size = 0.1,
          ylim = c(0, 50)) +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(size = 20, hjust = 1),
          axis.text.y = element_text(size = 20),
          panel.grid.major = element_line(size = 0.5, color = "gray90")) + 
    # ggsave(filename = file.path(outpath, experiment, str_c(bam_name, "IR_phase.png", sep = ".")), width = 10, height = 10)
    ggsave(filename = file.path(outpath, "test2.png"), width = 10, height = 10)
  
  # return
  return(bam_name)
  
}))




