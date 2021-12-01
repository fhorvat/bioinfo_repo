### INFO: 
### DATE: Tue Feb 04 16:44:06 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Roovers_2015_CellRep_GSE64942/Analysis/expression/compare")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# miRBase path
mirbase_path <- file.path(inpath, "../miRBase/miRBase.22.Btau_5.0.1.20200202.genBankseqnames.FPM_mean.csv")

# Ensembl annotated miRNAs
ensembl_path <- file.path(inpath, "../miRNA/ensembl.96.ARS-UCD1.2.20190809.UCSCseqnames.miRNA.FPM_mean.csv")
  
######################################################## READ DATA
# read miRBase FPMs
mirbase_tb <- readr::read_csv(mirbase_path)

# read ensembl FPMs
ensembl_tb <- readr::read_csv(ensembl_path)

######################################################## MAIN CODE
# clean and join ensembl 
both_tb <- 
  ensembl_tb %>% 
  dplyr::mutate(gene_name = str_replace(gene_name, "-mir-", "-miR-")) %>% 
  dplyr::select(gene_id = gene_name, coordinates, GV_ensembl = GV) %>% 
  dplyr::inner_join(., mirbase_tb %>% dplyr::select(gene_id, GV_mirbase = GV), by = "gene_id") 

# build linear regression model
linearMod <- lm(GV_mirbase ~ GV_ensembl, data = both_tb)  # build linear regression model on full data
summary(linearMod)

# join with ERCC fpkm, calculate log2 of FPKM and molecule numbers, prepare for plotting
plot_tb <- 
  both_tb %>% 
  dplyr::mutate(log2_mirbase = log2(GV_mirbase + 0.1), 
                log2_ensembl = log2(GV_ensembl + 0.1))

# plot without gene labels
ln_plot <- 
  ggplot() +
  geom_point(data = plot_tb, 
             mapping = aes(x = GV_mirbase, y = GV_ensembl), alpha = 1, size = 3) +
  geom_smooth(data = plot_tb,
              mapping = aes(x = GV_mirbase, y = GV_mirbase),
              method = lm) +
  scale_x_continuous(limits = c(0, 500)) +
  scale_y_continuous(limits = c(0, 500)) +
  xlab("log2(miRBase FPM)") +
  ylab("log2(ensembl FPM)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none")

# save
ggsave(filename = file.path(outpath, "linear_regression.compare_miRBase_and_ensembl.FPM.png"), plot = ln_plot, height = 10, width = 10)


