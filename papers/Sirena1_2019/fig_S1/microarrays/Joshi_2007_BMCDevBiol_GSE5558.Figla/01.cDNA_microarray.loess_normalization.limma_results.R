### INFO: 
### DATE: Fri Aug 16 17:50:40 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/diffExp/Joshi_2007_BMCDevBiol_GSE5558")

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

library(limma)
library(openxlsx)
library(ggrepel)
library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Joshi_2007_BMCDevBiol_GSE5558"

# set outpath
outpath <- getwd()

# annotation path
annotation_path <- file.path(inpath, "GPL3771.annot.tidy.csv")

# sample table path
sample_table_path <- file.path(inpath, "Joshi_2007_BMCDevBiol_GSE5558.sampleTable.csv")

# intensity values path
intensity_path <- list.files(inpath, "GSM[0-9]+\\.csv", full.names = T)

######################################################## READ DATA
# read annotation table
annotation_tb <- readr::read_csv(annotation_path)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
### prepare data for analysis
# filter table to include only newborn
sample_table_tidy <- 
  sample_table %>%
  dplyr::filter(stage == "newborn")

# get only samples from newborn 
intensity_filtered <- 
  sample_table_tidy %>% 
  dplyr::select(Label = geo_accession) %>% 
  dplyr::mutate(FileName = str_c(Label, ".csv")) %>% 
  as.data.frame(.)


### read intensity values
# weight function - if spot has Flag < 0 it get 0 weight
f <- function(x) as.numeric(x$Flags == 0)

# read RGList
rglist <- limma::read.maimages(files = intensity_filtered, 
                               path = inpath, 
                               names = intensity_filtered$Label, 
                               columns = c("R" = "CH1_MEAN", "G" = "CH2_MEAN", "Rb" = "CH1_BKD_Mean", "Gb" = "CH2_BKD_Mean"),
                               other.columns = c("Flags" = "Flags"), 
                               annotation = "ID_REF", 
                               sep = ",", 
                               wt.fun = f)

# correct background
rglist_corrected <- backgroundCorrect(rglist, method = "normexp", offset = 0)
?limma::backgroundCorrect
# normalize
rglist_normalized <- normalizeWithinArrays(rglist_corrected, method = "loess")
?limma::normalizeWithinArrays

# Estimate the fold changes and standard errors by fitting a linear model for each gene. 
# The design matrix indicates which arrays are dye-swaps
fit <- lmFit(rglist_normalized, design = c(-1, -1, -1, -1, 1, 1, 1, 1))

# Apply empirical Bayes smoothing to the standard errors.
fit <- eBayes(fit)

# get the results, tidy
results_tb <- 
  topTable(fit, n = nrow(rglist)) %>% 
  as_tibble(.) %>% 
  dplyr::select(ID = ID_REF, lfc = logFC, mean = AveExpr, padj = adj.P.Val) %>% 
  dplyr::left_join(., annotation_tb, by = "ID") %>% 
  dplyr::arrange(padj)


### save results
# save all results to xlsx file
results_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = results_wb, sheetName = "Figla_Null_vs_WT.all")
openxlsx::writeData(wb = results_wb, sheet = "Figla_Null_vs_WT.all", x = results_tb)
openxlsx::saveWorkbook(wb = results_wb, 
                       file = file.path(outpath, str_c("diffExp", "VMSRMus20K", "Figla_Null_vs_WT",
                                                       "newborn_ovaries", "limma",
                                                       "results.xlsx", sep = ".")), 
                       overwrite = TRUE)


### plot MA plot
# data for plot
plot_df <- 
  results_tb %>% 
  dplyr::select(mean, lfc, padj, gene_symbol, ID, gene_name, coordinates) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, padj > 0.05, "no")) %>%
  dplyr::mutate(regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
  dplyr::arrange(regulation)

# get annotation data
annotation_df <-
  plot_df %>%
  dplyr::filter(gene_symbol %in% c("C86187", "Nobox", "Lhx8", "Figla"))

# plot
ma_plot <- 
  ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), shape = 20, size = 5) +
  geom_point(data = annotation_df, aes(x = mean, y = lfc), color = "black", alpha = 1, size = 5, shape = 20) +
  # scale_x_continuous(limits = c(-4, 4), breaks = c(-4:4)) +
  scale_y_continuous(limits = c(-4, 4), breaks = c(-12:12)) +
  scale_colour_manual(values = c(no = "gray50", down = "#1a75ff", up = "red3")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none")

# save plot
ggsave(filename = file.path(outpath, str_c("diffExp", "VMSRMus20K", "Figla_Null_vs_WT",
                                           "newborn_ovaries", "limma",
                                           "MA_plot.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)

# add lables to genes
ma_plot <- 
  ma_plot + 
  geom_label_repel(data = annotation_df, aes(x = mean, y = lfc, label = gene_symbol), 
                   fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save plot
ggsave(filename = file.path(outpath, str_c("diffExp", "VMSRMus20K", "Figla_Null_vs_WT",
                                           "newborn_ovaries", "limma",
                                           "MA_plot.labels.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)


# ### check MA plots 
# pdf(file.path(outpath, "ma_plots.01.raw_arrays.pdf"), width = 10, height = 10)
# for(i in 1:8){
#   limma::plotMA(object = rglist, array = i)
# }
# dev.off()
# 
# pdf(file.path(outpath, "ma_plots.02.corrected_arrays.pdf"), width = 10, height = 10)
# for(i in 1:8){
#   limma::plotMA(object = rglist_corrected, array = i)
# }
# dev.off()
# 
# pdf(file.path(outpath, "ma_plots.03.normalized_arrays.pdf"), width = 10, height = 10)
# for(i in 1:8){
#   limma::plotMA(object = rglist_normalized, array = i)
# }
# dev.off()