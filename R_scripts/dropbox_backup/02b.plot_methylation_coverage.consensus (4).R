### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/methylation/mapped/Nov_2020/FLI_consensus.IAPLTR4_I")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- inpath

# methylation paths
meth_path_list <- list.files(inpath, ".*_bismark_bt2_pe\\.bismark\\.cov\\.gz", full.names = T) 

# sequences path
seq_path <- "../../../bismark_index/FLI_consensus.IAPLTR4_I"
seq_path <- list.files(seq_path, ".*\\.fasta$", full.names = T)

######################################################## READ DATA
# read table coverage data
meth_tb <- purrr::map(meth_path_list, function(path){
  
  # read table
  readr::read_delim(path, delim = "\t", col_names = c("repName", "start", "end", "percentage_meth", "count_meth", "count_unmeth")) %>% 
    dplyr::mutate(sample_id = basename(path) %>% str_remove(., "_bismark_bt2_pe\\.bismark\\.cov\\.gz"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(percentage_unmeth = 100 - percentage_meth) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO")))

# read consensus sequence 
cons_seq <- Biostrings::readDNAStringSet(seq_path)

######################################################## MAIN CODE
### plot coverage for different repNames
purrr::map(unique(meth_tb$repName), function(rmsk){
  
  # get percentage of methylation
  methylation_tb <- 
    meth_tb %>% 
    dplyr::filter(repName == rmsk) %>% 
    dplyr::select(repName, sample_id, genotype, pos = start, percentage_meth, percentage_unmeth) %>% 
    tidyr::pivot_longer(., -c(repName, sample_id, genotype, pos), names_to = "meth", names_prefix = "percentage_", values_to = "count") %>% 
    dplyr::filter(meth == "meth")
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_bar(data = methylation_tb,
             mapping = aes(x = pos, y = count), fill = "black", width = 10, stat = "identity") +
    facet_grid(rows = vars(genotype)) +
    xlab("") +
    ylab("") +
    scale_x_continuous(limits = c(0, width(cons_seq)), 
                       # labels = c(0, methylation_tb$pos, width(cons_seq)), 
                       # breaks = c(0, methylation_tb$pos, width(cons_seq))
                       ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  # save plot
  ggsave(plot = baseplot, 
         filename = file.path(outpath, str_c("IAPLTR4_I.FLI_consensus.methylation", "png", sep = ".")), 
         width = 15, height = 7.5)
  
  # return
  return(rmsk)
  
})
