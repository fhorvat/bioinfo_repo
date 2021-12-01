### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/testis.RNA_seq")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### 13dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression")

# FPKM table path
fpkm_tb_path.13dpp <- list.files(expression_path, ".*\\.FPKM_mean\\.csv$", full.names = T)

# results path
results_path.13dpp <- list.files(expression_path, ".*\\.all_results\\.xlsx$", full.names = T, recursive = T)


### 8.5dpp
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.8.5dpp.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression")

# FPKM table path
fpkm_tb_path.8.5dpp <- list.files(expression_path, ".*\\.FPKM_mean\\.csv$", full.names = T)

# results path
results_path.8.5dpp <- list.files(expression_path, ".*\\.significant_results\\.xlsx$", full.names = T, recursive = T)

######################################################## READ DATA
### 8.5dpp 
# read fpkm table
fpkm_tb.8.5dpp <- readr::read_csv(fpkm_tb_path.8.5dpp)

# read results table 
results.8.5dpp <- openxlsx::read.xlsx(results_path.8.5dpp) %>% as_tibble(.)


### 13dpp 
# read fpkm table
fpkm_tb.13dpp <- readr::read_csv(fpkm_tb_path.13dpp)

# read results table 
results.13dpp <- openxlsx::read.xlsx(results_path.13dpp) %>% as_tibble(.)

######################################################## MAIN CODE
# set name
table_name <- "MesAur1.RNA_seq"

### plot MA plot - 13dpp with 8.5 dpp guys with logFC > 5 in red
# get log2FC 
log2fc_tb <- 
  left_join(fpkm_tb.8.5dpp, fpkm_tb.13dpp, by = "gene_id") %>% 
  set_colnames(., str_remove(colnames(.), "Mov10l_")) %>% 
  # dplyr::filter_at(.vars = vars(contains("WT")), .vars_predicate = any_vars(. > fpkm_cutoff)) %>% 
  dplyr::mutate(log2FC_8.5dpp = log2(KO_8.5dpp / WT_8.5dpp), 
                log2FC_13dpp = log2(KO_13dpp / WT_13dpp), 
                log2FC_21dpp = log2(KO_21dpp / WT_21dpp)) %>% 
  dplyr::select(gene_id, log2FC_8.5dpp, log2FC_13dpp, log2FC_21dpp) %>% 
  dplyr::mutate_at(.vars = vars(contains("log2FC")), ~(replace(., is.nan(.), 0)))

# get genes with > 5 log2FC in 8.5 dpp
log2fc_8.5dpp_filt <- 
  log2fc_tb %>% 
  dplyr::filter(log2FC_8.5dpp > 6) 

# get table for MA plot
plot_df <-
  results.13dpp %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>%
  dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                padj = replace(padj, padj == 0, .Machine$double.xmin)) %>%
  dplyr::mutate(regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, padj > 0.1, "no"),
                regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
  dplyr::mutate(label = ifelse(gene_id %in% log2fc_8.5dpp_filt$gene_id, "yes", "no"),
                label = factor(label, levels = c("no", "yes"))) %>% 
  dplyr::arrange(label) 

# get axis limits
results_limits <-
  plot_df %>% 
  dplyr::summarise(x_limit = mean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = lfc %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.)) %>% 
  dplyr::mutate(y_limit = 15)

# plot
ma_plot <-
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = label, alpha = label), size = 5, shape = 20) +
  scale_x_log10(limits = c(0.01, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                     breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
  scale_colour_manual(labels = c(no = "log2FC 8.5dpp < 5", yes = "log2FC 8.5dpp > 5"),
                      values = c(no = "gray50", yes = "red2")) +
  scale_alpha_manual(values = c(no = 0.5, yes = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# save plot
ggsave(filename = file.path(outpath, 
                            str_c(table_name, "MA_plot.13dpp.with_8.5dpp_genes_labeled", "jpg", sep = ".")),
       plot = ma_plot, width = 10, height = 10)



