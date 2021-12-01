#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
### PATH: /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/SplicingAnalysis/PBS.SplicingAnalysis.TEMP.R
rm(list=ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/")

################################################################################### LIBRARIES
lib.path <- "/common/WORK/fhorvat/R_library/vfranke"
source(file.path(lib.path, "FileLoader.R"))
source(file.path(lib.path, "FormatConverters.R"))
source(file.path(lib.path, "BamWorkers.R"))
source(file.path(lib.path, "ScanLib.R"))
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)
library(genomation)
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)

################################################################################### PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/RData/Splicing_Analysis_RData"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/results"

tot_counts_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/documentation/CNOT6L_library_size.txt"

################################################################################### SCRIPT PARAMS
# register workers for parallel computation
registerDoMC(21)

################################################################################### READING TABLES
# total counts
sample_table <- read_delim(file = tot_counts_path, delim = "\t")
tot_counts <- as.integer(sample_table$library_size) %>% 
  set_names(., sample_table$name)

# normalize library size
tot_norm <- round(1e7 / tot_counts, 2)

# order of samples
ord <- c("s_GV_WT_X4", "s_GV_WT_X5", "s_GV_WT_X6", "s_GV_KO_X1", "s_GV_KO_X2", "s_GV_KO_X3",
         "s_MII_WT_X10", "s_MII_WT_X11",  "s_MII_WT_X12", "s_MII_KO_X7", "s_MII_KO_X8", "s_MII_KO_X9", 
         "s_1C_WT_X16", "s_1C_WT_X17", "s_1C_WT_X18", "s_1C_KO_X13", "s_1C_KO_X14", "s_1C_KO_X15", 
         "s_GV_Hamster_X19", "s_GV_Hamster_X20")

################################################################################### MAIN CODE
# read RData files
rfiles <- list.files(inpath, full.names = TRUE, pattern = "RData")
rfiles <- rfiles[!str_detect(rfiles, "X20|X19")]
l <- foreach(i = 1:length(rfiles))%dopar%{
  print(i)
  Assigner(rfiles[i], "a")
  return(a)
}

# set names, order samples
names(l) <-  
  left_join(tibble(ID = str_replace(basename(rfiles), "_.*", "")),
            sample_table, 
            by = "ID") %$% 
  name

ord <- ord[ord %in% names(l)]
l <- l[ord]
sample_names <- names(l)

# extract tables from list, add sample column to each data.frame in list
sl <- lapply(l, function(x) x[[2]])
rm(l); gc()

sl <- lapply(names(sl), function(x){
  a <- sl[[x]]
  a$samp <- x
  a
})

# bind data.frames in list to one data.frame, change sample order, change column names
sl <- rbindlist(sl)
sl$samp <- factor(sl$samp, levels = sample_names)
setnames(sl, 3:8, str_replace(colnames(sl)[3:8], "genes_", ""))
sl$ex[sl$int] <- FALSE

################################################################################### tables for plots and plots
### reads which are not expressed in exons/introns in MII stage and which are expressed in 1-cell stage introns
# sum counts of reads mapped to exon/intron for each sample
d <- sl[sl$sel_mii, list(ex = sum(ex), int = sum(int)), by = samp]
m <- 
  reshape2::melt(d, id.vars = "samp") %>% 
  dplyr::mutate(norm = value * tot_norm[samp]) %>% 
  as.data.table(.)

# plot count sums
ggplot(m, aes(x = samp, y = norm, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ggtitle("MII removed") +
  ggsave(paste0(outpath, "/CNOT6L_exonIntron_counts.pdf"), width = 9.93, height = 5.86)

# get mean count sums 
m_mean <- 
  m %>% 
  mutate(samp = str_replace(samp, "_X.*", "")) %>% 
  dplyr::group_by(samp, variable) %>% 
  dplyr::summarise(norm = mean(norm)) %>% 
  ungroup() %>% 
  dplyr::mutate(samp = factor(samp, levels = c("s_GV_WT",  "s_GV_KO",  "s_MII_WT", "s_MII_KO", "s_1C_WT", "s_1C_KO")))

# plot mean count sums
ggplot(m_mean, aes(x = samp, y = norm, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ggtitle("MII removed") +
  ggsave(paste0(outpath, "/CNOT6L_mean_exonIntron_counts.pdf"), width = 9.93, height = 5.86)

################################################################################### ratio
# calculate ratio of read counts mapped to introns / exons
m_rat <- 
  m[, norm[variable == "int"] / norm[variable == "ex"], by = samp] %>% 
  data.table::setnames(., 2, "ratio") %>% 
  dplyr::mutate(ratio = replace(ratio, is.infinite(ratio), 0))

# plot ratio
ggplot(m_rat, aes(x = samp, y = ratio)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ggtitle("MII removed, Intron/Exon ratio") +
  ggsave(paste0(outpath, "/CNOT6L_exonIntron_ratio.pdf"), width = 9.93, height = 5.86)

# calculate mean ratio
m_rat_mean <- 
  m_rat %>% 
  dplyr::mutate(samp = str_replace(samp, "_X.*", "")) %>% 
  dplyr::group_by(samp) %>% 
  dplyr::summarise(ratio = mean(ratio)) %>% 
  ungroup() %>% 
  dplyr::mutate(samp = factor(samp, levels = c("s_GV_WT",  "s_GV_KO",  "s_MII_WT", "s_MII_KO", "s_1C_WT", "s_1C_KO")))

# plot mean ratios
ggplot(m_rat_mean, aes(x = samp, y = ratio)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ggtitle("MII removed, Intron/Exon ratio") +
  ggsave(paste0(outpath, "/CNOT6L_mean_exonIntron_ratio.pdf"), width = 9.93, height = 5.86)

################################################################################### black and white ratio plot
# ggplot(m.rat, aes(x = samp, y = ratio)) +
#   geom_bar(stat = "identity", color = "black", fill = "black") +
#   scale_x_discrete(drop = FALSE) + 
#   theme_bw() + 
#   theme(panel.border = element_rect(size = .1, colour = "black"),
#         axis.ticks = element_line(size = 1),
#         axis.text = element_text(size = 20),
#         axis.text.x = element_text(angle = 90),
#         axis.title = element_text(size = 20),
#         panel.grid = element_blank(),
#         legend.key.size = unit(1, "cm"),
#         legend.text = element_text(size = 15),
#         legend.title = element_text(size = 15),
#         legend.position = "none") +
#   ylab("Unspliced / Spliced") +
#   xlab("") +
#   ggsave(paste0(outpath, "/CNOT6L_mean_exonIntron_ratio_black.pdf"), width = 9.93, height = 5.86)
