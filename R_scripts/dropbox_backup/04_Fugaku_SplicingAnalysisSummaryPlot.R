#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
### PATH: /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/SplicingAnalysis/PBS.SplicingAnalysis.TEMP.R
rm(list=ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/")

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

################################################################################### PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/Splicing_Analysis_RData"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/results"

tot_counts_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/spliced_reads/Fugaku/output/documentation/Fugaku_library_size.txt"

################################################################################### SCRIPT PARAMS
# register workers for parallel computation
registerDoMC(21)

################################################################################### READING TABLES
# total counts
tot_counts <- read_delim(file = tot_counts_path, delim = "\t")
tot_counts <- as.integer(tot_counts$library_size) %>% 
  set_names(., tot_counts$ID)

# normalize library size
tot_norm <- round(1e7 / tot_counts, 2)

# order of samples
ord <- c("s_GV.WE", "s_MII.PA", "s_MII.WE", 
         "s_1cell.PA", "s_1cell.WE", "s_1cell.WE_DNAm", 
         "s_2cell.WE", "s_2cell.WE_DNAm", 
         "s_4cell.WE", "s_Molura.WE", "s_Blast.WE")

################################################################################### MAIN CODE
# read RData files
rfiles <- list.files(inpath, full.names = TRUE, pattern = "RData")
l <- foreach(i = 1:length(rfiles))%dopar%{
  print(i)
  Assigner(rfiles[i], "a")
  return(a)
}

# set names, order samples
names(l) <- str_replace(basename(rfiles), ".SpliceDonAcc.RData", "")
ord <- ord[ord %in% names(l)]
l <- l[ord]

# extract tables from list, add sample column to each data.frame in list
sl <- lapply(l, function(x) x[[2]])
sl <- lapply(names(sl), function(x){
  a <- sl[[x]]
  a$samp <- x
  a
})

# bind data.frames in list to one data.frame, change sample order, change column names
sl <- rbindlist(sl)
sl$samp <- factor(sl$samp, levels = names(l))
setnames(sl, 3:8, str_replace(colnames(sl)[3:8], "genes_", ""))
sl$ex[sl$int] <- FALSE

################################################################################### tables for plots and plots
### reads which are not expressed in exons/introns in MII stage and which are expressed in 1-cell stage introns

# sum counts of reads mapped to exon/intron for each sample
d <- sl[sl$sel_mii, list(ex = sum(ex), int = sum(int)), by = samp]
m <- melt(d[!str_detect(d$samp, "PA"), ])
m$norm <- m$value * tot_norm[m$samp]

ggplot(m, aes(x = samp, y = norm, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ggtitle("MII removed") +
  ggsave(paste0(outpath, "/Fugaku_exonIntron_sumCounts.pdf"), width = 9.93, height = 5.86)


# plot ratio of read counts mapped to introns / exons
m.rat <- m[, norm[variable == "int"] / norm[variable == "ex"], by = samp]
setnames(m.rat, 2, "ratio")

ggplot(m.rat, aes(x = samp, y = ratio)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.5, size = 10)) + 
  ggtitle("MII removed, Intron/Exon ratio") +
  ggsave(paste0(outpath, "/Fugaku_exonIntron_ratio.pdf"), width = 9.93, height = 5.86)

ggplot(m.rat, aes(x = samp, y = ratio)) +
  geom_bar(stat = "identity", color = "black", fill = "black") +
  theme_bw() + 
  theme(panel.border = element_rect(size = .1, colour = "black"),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 20),
        panel.grid = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "none") +
  ylab("Unspliced / Spliced") +
  xlab("") +
  ggsave(paste0(outpath, "/Fugaku_exonIntron_ratio_black.pdf"), width = 9.93, height = 5.86)
