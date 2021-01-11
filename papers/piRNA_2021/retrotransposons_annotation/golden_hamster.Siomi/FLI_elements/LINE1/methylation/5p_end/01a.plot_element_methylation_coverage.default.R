### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/methylation/mapped/Nov_2020/5p_end")

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
library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- inpath

# bigWigs path
bw_path <- list.files(inpath, ".*bw$", full.names = T) 
bw_path <- bw_path[!str_detect(bw_path, "PE_s|raw_coverage")]

# coverage paths
coverage_path <- list.files(inpath, ".*_bismark_bt2_pe\\.filt\\.bismark\\.cov\\.gz", full.names = T) 

# sequences path
seq_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/LINE/methylation/bismark_index/5p_sequences"
seq_path <- list.files(seq_path, ".*\\.fasta$", full.names = T)

######################################################## READ DATA
# read table coverage data
meth_tb <- purrr::map(coverage_path, function(path){
  
  # read table
  readr::read_delim(path, delim = "\t", col_names = c("repName", "start", "end", "percentage_meth", "count_meth", "count_unmeth")) %>% 
    dplyr::mutate(sample_id = basename(path) %>% str_remove(., "_bismark_bt2_pe.filt.bismark.cov.gz"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(percentage_unmeth = 100 - percentage_meth) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO")))

# read fasta 
fasta_seq <- Biostrings::readDNAStringSet(seq_path)

######################################################## MAIN CODE
### plot coverage for different repNames
purrr::map(unique(meth_tb$repName), function(rmsk){
  
  # get percentage of methylation
  methylation_tb <- 
    meth_tb %>% 
    dplyr::filter(repName == rmsk) %>% 
    dplyr::select(repName, sample_id, genotype, pos = start, percentage_meth, percentage_unmeth) %>% 
    tidyr::pivot_longer(., -c(repName, sample_id, genotype, pos), names_to = "meth", names_prefix = "percentage_", values_to = "count") %>% 
    dplyr::mutate(meth = factor(meth, levels = c("unmeth", "meth")))
  
  # get total number of reads
  nreads_tb <- 
    meth_tb %>% 
    dplyr::filter(repName == rmsk) %>% 
    dplyr::mutate(count_total = count_meth + count_unmeth) %>% 
    dplyr::select(repName, sample_id, genotype, pos = start, count_total)
  
  # get sequence
  rmsk_seq <- 
    fasta_seq[[rmsk]] %>% 
    as.character(.)
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_bar(data = methylation_tb,
             mapping = aes(x = pos, y = count, fill = meth), width = 0.8, stat = "identity") +
    geom_text(data = nreads_tb, 
              mapping = aes(x = pos, y = 50, label = count_total), size = 3) + 
    facet_grid(rows = vars(genotype)) +
    xlab("") +
    ylab("") +
    scale_x_continuous(limits = c(0, nchar(rmsk_seq)), breaks = 1:nchar(rmsk_seq), labels = unlist(str_split(rmsk_seq, ""))) +
    # scale_y_continuous(limits = c(0, 90)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.title = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(plot = baseplot, 
         filename = file.path(outpath, str_c("methylation_coverage", rmsk, "bw", "png", sep = ".")), 
         width = 15, height = 7.5)
  
  # return
  return(rmsk)
  
})
