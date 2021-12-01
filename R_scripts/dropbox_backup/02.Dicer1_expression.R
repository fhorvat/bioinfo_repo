### INFO: 
### DATE: Thu Apr 08 19:22:10 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/tissues.Dicer_MT_HET.2021_Apr/Analysis/Dicer_expression")

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

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/rat/rn6.Rnor_6.0.GCA_000001895.4"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# tissue expression path
expression_path_1 <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/tissues.Dicer_MT_HET.2021_Apr/Analysis/expression"
expression_path_1 <- file.path(expression_path_1, "ensembl.99.Rnor_6.0.20200415.UCSCseqnames.FPKM.csv")

# oocyte expression path
expression_path_2 <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_rat/datasets/oocytes.Dicer_MT_HET.2021_Apr/Analysis/expression"
expression_path_2 <- file.path(expression_path_2, "ensembl.99.Rnor_6.0.20200415.UCSCseqnames.FPKM.csv")

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read expression tables
expression_tb_1 <- readr::read_csv(expression_path_1)
expression_tb_2 <- readr::read_csv(expression_path_2)

######################################################## MAIN CODE
# join tables
expression_tb <- 
  expression_tb_1 %>% 
  dplyr::left_join(., expression_tb_2, by = "gene_id") %>% 
  dplyr::select(gene_id, starts_with("s_")) %>% 
  tidyr::pivot_longer(data = ., cols = -gene_id, names_to = "sample_id", values_to = "fpkm") %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "MT_HET|WT"), 
                tissue = str_remove_all(sample_id, "^s_|_r1.SE") %>% 
                  str_remove(., genotype) %>% 
                  str_remove(., "_") %>% 
                  str_replace(., "GVDicer_", "oocyte")) %>% 
  dplyr::left_join(., genes_info, by = "gene_id") %>% 
  dplyr::filter(gene_name == "Dicer1") 

# save as table
dicer_tb <- 
  expression_tb %>% 
  dplyr::select(tissue, genotype, fpkm) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "MT_HET"))) %>% 
  dplyr::arrange(genotype) %>% 
  tidyr::pivot_wider(id_cols = tissue, names_from = genotype, values_from = fpkm) %>% 
  readr::write_csv(., file = file.path(outpath, "ensembl.99.Rnor_6.0.20200415.Dicer1.FPKM.csv"))

# create plot table
plot_tb <- 
  dicer_tb %>% 
  tidyr::pivot_longer(cols = -tissue, names_to = "genotype", values_to = "fpkm") %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("WT", "MT_HET")))
  
# plot as barplot
dicer_barplot <- 
  ggplot() + 
  geom_bar(data = plot_tb, 
           mapping = aes(x = tissue, y = fpkm, fill = genotype), width = 0.8, stat = "identity", position = "dodge") + 
  # scale_fill_manual(values = c(miRNA = "#70ad47", mRNA = "#ffc000", repeats = "#ff0000", other_mapped = "#000000")) + 
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# save
ggsave(plot = dicer_barplot, filename = file.path(outpath, "barplot.ensembl.99.Rnor_6.0.20200415.Dicer1.FPKM.png"), width = 10, height = 10)


  
    
    
    
    
    
