library("dplyr")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("geneplotter")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("scales")
library(readr)
# show_col()

options(bitmapType = 'cairo')
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Analysis")

################################################################## functions
# calculate and plot PCA
PCAplot <- function(data_df, sample_table_df, plot_name){
  
  #### PCA
  # filteres and orders data by highest variance in first 500 rows 
  ntop <- 500
  rv <- genefilter::rowVars(data_df) 
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  # calculates pca
  pca <- prcomp(t(data_df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  # makes data.frame
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  name = sample_table_df$"name", 
                  Treatment.Control = sample_table_df$"Treatment.Control", 
                  lib_size = sample_table_df$"lib_size", 
                  Time.Course = sample_table_df$"Time.Course")
  d$Treatment.Control <- factor(d$Treatment.Control, levels = c("WT", "KO", "Lnc5", "Lnc12", "Lnc21"))
  
  # plot
  # show_col(hue_pal()(8))
  ggplot(d, aes(PC1, PC2)) + 
    geom_point(aes(size = lib_size, color = Treatment.Control)) +
    geom_text(aes(label = name), vjust = "inward", hjust = "inward", nudge_y = -0.2, nudge_x = 0.2, size = 2, color = "black") +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    scale_color_manual(values = color_pallete) + 
    ggsave(plot_name)
  
  return(plot_name)
}

################################################################## reading data
# lnc_2016 experiment sample table
sample_table_lnc2016 <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/RNAseq_2016_11_23_sampleTable.csv", header = T) %>%
  dplyr::select(ID = sample, Time.Course = Sample, Treatment.Control = Type) %>%
  mutate(Treatment.Control = gsub(" B6", "", Treatment.Control), 
         ID = gsub("_16.*|_[A,C,T,G].*", "", ID), 
         name = paste(gsub(".*_[1-8]_", "", ID), Treatment.Control, sep = "_"), 
         experiment = "lnc_Nov2016") %>% 
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T,
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_16.*|_[A,C,T,G].*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2]) / 10^6)), 
            by = "ID") %>% 
  mutate(name = ifelse(grepl("HV", ID), paste0(name, "_new"), name))

# CNOT6L sample table
sample_table_CNOT6L <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_sample_list_11919R_2015_10_29.csv", header = T) %>%
  dplyr::select(ID, Time.Course, Treatment.Control) %>%
  mutate(name = paste(ID, Time.Course, Treatment.Control, sep = "_"),  
         experiment = "CNOT6L") %>%
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T, 
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_.*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2]) / 10^6)), 
            by = "ID")

# combine_tables
sample_table_all <- rbind(sample_table_lnc2016, sample_table_CNOT6L)
sample_table_all <- sample_table_all[-c(30:31), ]

################################################################## reading FPKMs
fpkm_all <- 
  read_csv("fpkm_knownGene_mm10_lnc2016_CNOT6L.csv") %>% 
  dplyr::select(-matches("Hamster|entrez|width")) %>% 
  as.data.frame(.)

################################################################## PCA and plot
# check if fpkm columns are parallel with sample data.frame rows
all(colnames(fpkm_all) == sample_table_all$name)

# mannualy set colors to the stages
color_pallete <- hue_pal()(8)[c(1, 5, 7, 6, 3)]
names(color_pallete) <- c("WT", "KO", "Lnc5", "Lnc12", "Lnc21")

# # old plots
# - CNOT6L all
# - CNOT6L GV (WT + KO)
# - lnc2016 old
# - lnc2016 old + CNOT6L all
# - lnc2016 old + CNOT6L GV (WT + KO)
# - lnc2016 old + CNOT6L GV WT
# 
# # new plots
# - lnc2016 old and new all c(1:11)
# - lnc2016 old and new lnc21 + WT 
# - lnc2016 old lnc5 + lnc12, new lnc21 + WT
# 
# - lnc2016 old and new + CNOT6L all
# - lnc2016 old and new + CNOT6L GV (WT + KO)
# - lnc2016 old and new + CNOT6L GV WT
# 
# - lnc2016 new + CNOT6L all
# - lnc2016 new + CNOT6L GV (WT + KO)
# - lnc2016 new + CNOT6L GV WT

sample_filter <- c(1:11)
data_df <- fpkm_all[, sample_filter]
sample_table_df <-  sample_table_lnc2016[sample_filter, ]
plot_name <- "test.png"
PCAplot(data_df, sample_table_df, plot_name)
