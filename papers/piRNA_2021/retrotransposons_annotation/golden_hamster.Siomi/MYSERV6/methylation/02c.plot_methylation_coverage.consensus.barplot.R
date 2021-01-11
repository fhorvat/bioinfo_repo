### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/MYSERV6/methylation/mapped/Nov_2020/FLI_consensus")

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
seq_path <- "../../../bismark_index/FLI_consensus"
seq_path <- list.files(seq_path, ".*\\.fasta$", full.names = T)

######################################################## READ DATA
# read table coverage data
meth_tb <- purrr::map(meth_path_list, function(path){
  
  # read table
  readr::read_delim(path, delim = "\t", col_names = c("repName", "start", "end", "percentage_meth", "count_meth", "count_unmeth")) %>% 
    dplyr::mutate(sample_id = basename(path) %>% str_remove(., "_bismark_bt2_pe\\.bismark\\.cov\\.gz"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_HET|Mov10l1_KO")) %>% 
  dplyr::filter(genotype != "Mov10l1_WT") %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Mov10l1_HET", "Mov10l1_KO"))) %>% 
  dplyr::filter(count_meth + count_unmeth >= 4)

# read sequence
fli_seq <- Biostrings::readDNAStringSet(seq_path)

######################################################## MAIN CODE
# get percentage of methylation
methylation_tb <- 
  meth_tb %>% 
  dplyr::select(genotype, pos = start, percentage_meth) %>% 
  tidyr::pivot_wider(id_cols = pos, values_from = percentage_meth, names_from = genotype)  %>% 
  dplyr::filter(!is.na(Mov10l1_HET), !is.na(Mov10l1_KO)) %>% 
  dplyr::slice_min(pos, n = 20) %>% 
  dplyr::mutate(rank = order(pos)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_longer(cols = -c(pos, rank), names_to = "genotype", values_to = "meth") %>% 
  dplyr::mutate(unmeth = 100 - meth) %>% 
  tidyr::pivot_longer(cols = meth:unmeth, names_to = "meth", values_to = "count") %>% 
  dplyr::mutate(meth = factor(meth, levels = c("unmeth", "meth"))) %>% 
  dplyr::mutate(rank = str_c("r_", rank), 
                rank = factor(rank, levels = str_c("r_", 1:20))) %>% 
  arrange(rank) %>% 
  dplyr::mutate(rank = factor(str_c(rank, ".pos_", pos), levels = unique(str_c(rank, ".pos_", pos))))

# create plot
baseplot <- 
  ggplot() +
  geom_bar(data = methylation_tb,
           mapping = aes(x = rank, y = count, fill = meth), width = 0.8, stat = "identity") +
  scale_fill_manual(values = c(meth = "black", unmeth = "gray80")) +
  facet_grid(rows = vars(genotype)) +
  xlab("") +
  ylab("") +
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
       filename = file.path(outpath, str_c("MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.top_10.consensus.methylation.barplot", "4_reads_filt", "png", sep = ".")), 
       width = 15, height = 7.5)


### create position plot
# create plot
baseplot <- 
  ggplot() +
  geom_bar(data = methylation_tb,
           mapping = aes(x = pos, y = count, fill = meth), width = 0.8, stat = "identity") +
  scale_fill_manual(values = c(meth = "black", unmeth = "gray80")) +
  # scale_x_continuous(limits = c(1, max(methylation_tb$pos) + 10), 
  #                    breaks = seq(100, max(methylation_tb$pos) + 100, 100), 
  #                    labels = seq(100, max(methylation_tb$pos) + 100, 100)) + 
  scale_x_continuous(limits = c(1, nchar(fli_seq)), 
                     breaks = seq(100, nchar(fli_seq), 1000), 
                     labels = seq(100, nchar(fli_seq), 1000)) + 
  facet_grid(rows = vars(genotype)) +
  xlab("") +
  ylab("") +
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
       filename = file.path(outpath, str_c("MYSERV6-int.longer_than_4kb.RPM_testis_8.5dpp.ORFs.top_10.consensus.methylation.position_barplot", "4_reads_filt", "png", sep = ".")), 
       width = 15, height = 7.5)


