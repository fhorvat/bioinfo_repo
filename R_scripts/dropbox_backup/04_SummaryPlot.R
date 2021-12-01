#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 6. 4. 2017.  
### AUTHOR: Filip Horvat
### PATH: 
rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads")

################################################################################### LIBRARIES
################################################################################### 
library(doMC)
library(GenomicAlignments)
library(genomation)
library(data.table)
library(stringr)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)

################################################################################### PATH VARIABLES
################################################################################### 
lib_path <- "/common/WORK/fhorvat/R_library/vfranke"

inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/RData/03_SplicingAnalysis"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/CNOT6L_Fugaku_together/results"

tot_counts_path_CNOT6L <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/documentation/CNOT6L_library_size.txt"
tot_counts_path_Fugaku <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/documentation/Fugaku_library_size.txt"

################################################################################### SOURCE FILES
################################################################################### 
source(file.path(lib_path, "FileLoader.R"))
source(file.path(lib_path, "FormatConverters.R"))
source(file.path(lib_path, "BamWorkers.R"))
source(file.path(lib_path, "ScanLib.R"))

################################################################################### FUNCTIONS
################################################################################### 
# plot 
plotIntronExon <- function(mii_transcripts_only){
  
  # get and summarize data
  if(mii_transcripts_only){
    
    # sum counts of reads mapped to exon/intron for each sample
    d <- sl[sl$genes_sel_mii, list(ex = sum(ex), int = sum(int)), by = samp]
    m <- 
      reshape2::melt(d, id.vars = "samp") %>% 
      dplyr::mutate(norm = value * tot_norm[samp]) %>% 
      as.data.table(.)
    plot_name <- "MIITrans"
    
  }else{
    
    d <- sl[, list(ex = sum(ex), int = sum(int)), by = samp]
    m <- 
      reshape2::melt(d, id.vars = "samp") %>% 
      dplyr::mutate(norm = value * tot_norm[samp]) %>% 
      as.data.table(.)
    plot_name <- "allTrans"
    
  }
  
  # plot count sums
  ggplot(m, aes(x = samp, y = norm, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_x_discrete(drop = FALSE) + 
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
    ggsave(paste0(outpath, "/Fugaku_CNOT6L_", plot_name, "_exonIntron_counts.pdf"), width = 9.93, height = 5.86)
  
  # get mean count sums 
  m_mean <- 
    m %>% 
    mutate(samp = str_replace(samp, "_X.*", "_CNOT6L")) %>% 
    dplyr::group_by(samp, variable) %>% 
    dplyr::summarise(norm = mean(norm)) %>% 
    ungroup() %>% 
    dplyr::mutate(samp = factor(samp, levels = ord_mean[ord_mean %in% samp]))
  
  # plot mean count sums
  ggplot(m_mean, aes(x = samp, y = norm, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_x_discrete(drop = FALSE) + 
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
    ggsave(paste0(outpath, "/Fugaku_CNOT6L_", plot_name, "_mean_exonIntron_counts.pdf"), width = 9.93, height = 5.86)
  
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
    ggtitle("Intron/Exon ratio") +
    ggsave(paste0(outpath, "/Fugaku_CNOT6L_", plot_name, "_exonIntron_ratio.pdf"), width = 9.93, height = 5.86)
  
  # calculate mean ratio
  m_rat_mean <- 
    m_rat %>% 
    dplyr::mutate(samp = str_replace(samp, "_X.*", "_CNOT6L")) %>% 
    dplyr::group_by(samp) %>% 
    dplyr::summarise(ratio = mean(ratio)) %>% 
    ungroup() %>% 
    dplyr::mutate(samp = factor(samp, levels = ord_mean[ord_mean %in% samp]))
  
  # plot mean ratios
  ggplot(m_rat_mean, aes(x = samp, y = ratio)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(drop = FALSE) + 
    theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
    ggtitle("Intron/Exon ratio") +
    ggsave(paste0(outpath, "/Fugaku_CNOT6L_", plot_name, "_mean_exonIntron_ratio.pdf"), width = 9.93, height = 5.86)
  
  cat(plot_name, "plot done \n")
  
}

################################################################################### SCRIPT PARAMS
################################################################################### 
# register workers for parallel computation
registerDoMC(21)

################################################################################### TABLES
################################################################################### 
# total counts
sample_table <- rbind(read_delim(file = tot_counts_path_Fugaku, delim = "\t") %>% 
                        dplyr::mutate(name = ID) %>% 
                        dplyr::select(ID, name, library_size), 
                      read_delim(file = tot_counts_path_CNOT6L, delim = "\t"))

tot_counts <- as.integer(sample_table$library_size) %>% 
  set_names(., sample_table$name)

# normalize library size
tot_norm <- round(1e7 / tot_counts, 2)

# order of samples
ord <- c("s_GV.WE", "s_GV_WT_X4", "s_GV_WT_X5", "s_GV_WT_X6", "s_GV_KO_X1", "s_GV_KO_X2", "s_GV_KO_X3",
         "s_MII.PA", "s_MII.WE", "s_MII_WT_X10", "s_MII_WT_X11",  "s_MII_WT_X12", "s_MII_KO_X7", "s_MII_KO_X8", "s_MII_KO_X9", 
         "s_1cell.PA", "s_1cell.WE", "s_1C_WT_X16", "s_1C_WT_X17", "s_1C_WT_X18", "s_1C_KO_X13", "s_1C_KO_X14", "s_1C_KO_X15", "s_1cell.WE_DNAm",
         "s_2cell.WE", "s_2cell.WE_DNAm", 
         "s_4cell.WE", 
         "s_Molura.WE", 
         "s_Blast.WE",
         "s_GV_Hamster_X19", "s_GV_Hamster_X20")

ord_mean <- c("s_GV.WE",  "s_GV_WT_CNOT6L",  "s_GV_KO_CNOT6L",
              "s_MII.PA", "s_MII.WE", "s_MII_WT_CNOT6L", "s_MII_KO_CNOT6L",
              "s_1cell.PA", "s_1cell.WE", "s_1C_WT_CNOT6L", "s_1C_KO_CNOT6L", "s_1cell.WE_DNAm", 
              "s_2cell.WE", "s_2cell.WE_DNAm", 
              "s_4cell.WE", 
              "s_Molura.WE", 
              "s_Blast.WE")

################################################################################### MAIN CODE
################################################################################### 
# read RData files
rfiles <- list.files(inpath, full.names = TRUE, pattern = "RData")
rfiles <- rfiles[!str_detect(rfiles, "X20|X19|PA")]

# get file names
rfiles_names <-
  tibble(rfiles = basename(rfiles)) %>%
  dplyr::mutate(ID = ifelse(str_detect(rfiles, "s_"),
                str_replace(rfiles, ".SpliceDonAcc.RData", ""), 
                str_replace(rfiles, "_.*", ""))) %>% 
  left_join(., sample_table, by = "ID") %$% 
  name

# load .RData files in parallel
l <- foreach(i = 1:length(rfiles))%dopar%{
  print(i)
  Assigner(rfiles[i], "a")
  return(a)
}

# set names, order samples
names(l) <- rfiles_names
ord <- ord[ord %in% names(l)]
l <- l[ord]
sample_names <- names(l)

# add sample column, bind do data.frame
sl <- lapply(sample_names, function(x){
  a <- l[[x]]
  a$samp <- x
  return(a)
}) %>% 
  rbindlist(.) %>% 
  dplyr::mutate(samp = factor(samp, levels = sample_names), 
                ex = replace(ex, int == T, F)) %>% 
  as.data.table(.)

rm(l); gc()

################################################################################### plot
plotIntronExon(mii_transcripts_only = T)
plotIntronExon(mii_transcripts_only = F)

