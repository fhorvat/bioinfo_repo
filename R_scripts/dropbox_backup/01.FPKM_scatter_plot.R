### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/ESC_vs_3T3_scatter")

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

library(openxlsx)
library(scales)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### genomic files
# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene GO info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.GOterms.csv$"), full.names = T)

  
### 3T3
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Darman_2015_CellRep_GSE72790"

# expression path
expression_path <- file.path(base_path, "Analysis/expression")

# FPKM table path
fpkm_tb_path.3T3 <- list.files(expression_path, ".*\\.FPKM_mean\\.csv$", full.names = T)


### mESC
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Freimer_2018_CurrBiol_GSE92761"

# expression path
expression_path <- file.path(base_path, "Analysis/expression")

# FPKM table path
fpkm_tb_path.mESC <- list.files(expression_path, ".*\\.FPKM_mean\\.csv$", full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

### 3T3 
# read fpkm table
fpkm_tb.3T3 <- 
  readr::read_csv(fpkm_tb_path.3T3) %>% 
  dplyr::select(gene_id, FPKM_3T3 = WT, gene_name, gene_biotype, gene_description)

### mESC 
# read fpkm table
fpkm_tb.mESC <- 
  readr::read_csv(fpkm_tb_path.mESC) %>% 
  dplyr::select(gene_id, FPKM_mESC = ESC)


######################################################## MAIN CODE
# set name
table_name <- "mESC_vs_3T3.RNA_seq"

# set FPKM cutoff
fpkm_cutoff <- 10

### plot crosshair plot - 3T3 vs. mESC
fpkm_tb_plot <- 
  dplyr::left_join(fpkm_tb.3T3, fpkm_tb.mESC, by = "gene_id") %>% 
  dplyr::left_join(genes_info, by = "gene_id") %>% 
  dplyr::select(gene_id, FPKM_3T3, FPKM_mESC, gene_name, gene_biotype, gene_description, go_id, go_name) %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  # dplyr::filter_at(.vars = vars(contains("FPKM")), .vars_predicate = any_vars(. > fpkm_cutoff)) %>% 
  dplyr::mutate(innate_immunity.GO = ifelse(str_detect(go_name, "interferon"), "interferon", "no_interferon"), 
                innate_immunity.GO = replace(innate_immunity.GO, is.na(innate_immunity.GO), "no_interferon"), 
                innate_immunity.GO = factor(innate_immunity.GO, levels = c("interferon", "no_interferon"))) %>% 
  dplyr::arrange(desc(innate_immunity.GO))

# scatter plot
scatter_plot <-
  ggplot(fpkm_tb_plot, aes(x = FPKM_3T3, y = FPKM_mESC, colour = innate_immunity.GO)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(values = c(no_interferon = "gray50", interferon = "red2")) +
  scale_alpha_manual(values = c(no_interferon = 0.5, interferon = 1)) +
  scale_x_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab(str_c("3T3 mean FPKM (Darman)")) +
  ylab(str_c("mESC mean FPKM (Freimer)")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c(table_name, "3T3_vs_mESC.FPKM.scatterplot.log", "png", sep = ".")),
       plot = scatter_plot,
       width = 10, height = 10)


