### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/MYSERV6/count_multimappers")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# bam subset path
experiment_name <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"
bam_paths <- file.path(inpath, "bam_subset", experiment_name)
bam_paths <- list.files(bam_paths, "s_testis_Mov10l1_KO_8.5dpp.*\\.bam$", full.names = T)

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.Siomi.20200701.clean.fa.out.gz")

######################################################## READ DATA
# # read clean repeatMasker
# rmsk_clean <- readr::read_delim(rmsk_clean_path, delim = "\t")
# 
# # subset MYSERV6-int
# subset_gr <- 
#   rmsk_clean %>% 
#   dplyr::filter(repName == "MYSERV6-int") %>% 
#   GRanges(.)

# read 
reads_tb <- purrr::map(bam_paths, function(path){
  
  # read alignments from bam
  chunk <- 
    readGAlignmentsList(path, 
                        param = ScanBamParam(what = "qname", 
                                             tag = c("NH")))
  
  # to table
  multi_tb <- 
    chunk %>% 
    unlist(.) %>% 
    as_tibble(.) %>% 
    dplyr::select(read_id = qname, count = NH) %>% 
    dplyr::distinct(.) %>% 
    dplyr::mutate(sample_id = path %>% basename(.) %>% str_remove(., "\\.bam"))
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
# calculate mean
reads_tb %>%
  # dplyr::group_by(sample_id) %>% 
  dplyr::summarize(multimap_mean = mean(count),
                   multimap_med = median(count))

# plot as histogram
multi_boxplot <- 
  ggplot() +
  geom_histogram(data = reads_tb, aes(x = count), binwidth = 10, colour = "black", fill = "white") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())

# save plot
ggsave(plot = multi_boxplot, 
       filename = file.path(outpath, str_c("MYSERV6-int", experiment_name, "multimappers", "histogram.png", sep = ".")), 
       width = 10, 
       height = 10)


# # plot as boxplots
# multi_boxplot <- 
#   ggplot() +
#   stat_boxplot(data = reads_tb, aes(x = sample_id, y = count), geom = "errorbar", size = 1.5) +
#   geom_boxplot(data = reads_tb, aes(x = sample_id, y = count), outlier.colour = NULL, outlier.shape = NA, size = 1.5) +
#   geom_jitter(data = reads_tb, aes(x = sample_id, y = count), alpha = 0.4, 
#               colour = "black", size = 1.5, width = 0.1, height = 0, show.legend = F) +
#   theme_bw(base_size = 10) +
#   theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         plot.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks.x = element_blank())
# 
# # save plot
# ggsave(plot = multi_boxplot, 
#        filename = file.path(outpath, str_c("MYSERV6-int", experiment_name, "multimappers", "boxplot.png", sep = ".")), 
#        width = 10, 
#        height = 10)





