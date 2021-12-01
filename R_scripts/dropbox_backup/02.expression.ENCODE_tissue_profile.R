### INFO: expression of lnc1 in ENCODE data set 
### DATE: Tue Aug 21 14:25:50 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/ovomucin_KO/Analysis/2021_paper/expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# path to ENCODE mapped samples
analysis_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417"

# sample table path
sample_table_path <- file.path(analysis_path, "Data/Documentation")
sample_table_path <- list.files(sample_table_path, ".*\\.sampleTable\\.csv", full.names = T)

# fpkm path
fpkm_path <- file.path(analysis_path, "Analysis/expression")
fpkm_path <- list.files(fpkm_path, ".*\\.FPKM_stats\\.csv", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
# gene ENSEMBL ID
gene_ensembl <- "ENSMUSG00000090891"

# set tissue order
tissue_order <- c("cns.E11.5", "cns.E14", "cns.E18", "frontallobe", "cortex", "cerebellum",
                  "stomach", "liver", "duodenum", "smintestine", "lgintestine", "colon", 
                  "lung", "heart", "bladder", "kidney", "thymus", "mammarygland", "spleen", 
                  "ovary", "testis", "placenta")

# tissue order with organ systems
tissue_order_df <- 
  tibble(tissue = tissue_order, 
         organ_system = c(rep("nervous", 6), 
                          rep("digestive", 6), 
                          rep("other", 7), 
                          rep("sex", 3)))

# get FPKM
fpkm_gene <- 
  fpkm_tb %>% 
  dplyr::filter(gene_id == gene_ensembl) %>% 
  dplyr::select(-c(coordinates, strand, gene_biotype, gene_description)) %>% 
  dplyr::left_join(., tissue_order_df, by = "tissue") %>% 
  dplyr::mutate(tissue = factor(tissue, levels = tissue_order)) %>% 
  dplyr::arrange(tissue)

# save
readr::write_csv(fpkm_gene, 
                 file.path(outpath, 
                           str_c("ENCODE_tissue_profile.barplot.FPKM_mean", gene_ensembl, "standard_deviation", "csv", sep = ".")))

# visualize
barplot_viridis <- 
  ggplot(fpkm_gene, 
         aes(x = tissue, y = avg_fpkm, fill = tissue)) + 
  geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin = avg_fpkm - SD, ymax = avg_fpkm + SD), width = 0.2) +
  geom_errorbar(aes(ymin = avg_fpkm - SE, ymax = avg_fpkm + SE), width = 0.2) +
  scale_fill_viridis_d() +
  ylab("FPKM") + 
  xlab("tissue") +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(barplot_viridis, 
       filename = file.path(outpath, 
                            str_c("ENCODE_tissue_profile.barplot.FPKM_mean", gene_ensembl, "standard_error", "png", sep = ".")), 
       width = 15, height = 10)


