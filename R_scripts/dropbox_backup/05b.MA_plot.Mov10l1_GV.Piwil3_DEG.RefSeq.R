### INFO: 
### DATE: Thu Jun 04 14:22:00 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Piwil1_KO_analysis")

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
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# ### genome files
# # gene info path
# genes_info_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed/annotation/Liftoff/MesAur1/RefSeq/hamster.sequel.draft-20200302.arrow.GCF_000349665.1_MesAur1.0.liftoff.geneInfo.csv"

# # read genes info
# genes_info <- readr::read_csv(genes_info_path)


### Mov10l1 KO path
# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq"

# expression path
expression_path <- file.path(base_path, "Analysis/expression.Siomi_RefSeq")

# results path
results_path.Mov10l1 <- list.files(expression_path, "all_biotype.*\\.all_results\\.xlsx$", full.names = T, recursive = T)


### Piwil3
# results path
results_path.Piwil3 <- file.path(inpath, "Piwil3 DEGs.xlsx")

######################################################## READ DATA
### Mov10l1 
# read results table 
results.Mov10l1 <- openxlsx::read.xlsx(results_path.Mov10l1, sheet = "Mov10l_KO_vs_Mov10l_WT") %>% as_tibble(.)

### Piwil
# read results table 
results.Piwil3 <- openxlsx::read.xlsx(results_path.Piwil3) %>% as_tibble(.)

######################################################## MAIN CODE
# clean table
results.Piwil3_tidy <- 
  results.Piwil3 %>% 
  dplyr::select(gene_name = Gene.Name, lfc_piwil3_KO = M.Value)

# join with our results
results.Piwil3_tidy %>% 
  dplyr::left_join(., results.Mov10l1 %>% dplyr::select(gene_name, lfc_mov10l1_KO = log2FoldChange, padj_mov10l1_KO = padj), by = "gene_name") %T>% 
  readr::write_csv(., file = file.path(outpath, "Piwil3_DEGs.Mov10l1_logFC_added.RefSeq.csv")) 

# set p-adjusted cut-off
padj_cutoff <- 0.01

# get significant result
results_df_sign <-
  results.Mov10l1 %>%
  dplyr::filter(padj <= padj_cutoff) %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::select(gene_name, regulation)


## MA plot
# data for plot
plot_df <-
  results.Mov10l1 %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name) %>%
  dplyr::left_join(., results.Piwil3_tidy, by = "gene_name") %>% 
  dplyr::left_join(., results_df_sign, by = "gene_name") %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                padj = replace(padj, padj == 0, .Machine$double.xmin)) %>% 
  dplyr::mutate(regulation = replace(regulation, is.na(regulation), "no"),
                regulation = replace(regulation, !is.na(lfc_piwil3_KO), "Piwil3"), 
                regulation = factor(regulation, levels = c("no", "up", "down", "Piwil3"))) %>%
  dplyr::arrange(regulation)

# get axis limits
results_limits <-
  plot_df %>% 
  dplyr::summarise(x_limit = mean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = lfc %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# plot
ma_plot <-
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 5, shape = 20) +
  scale_x_log10(limits = c(0.01, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                     breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated", Piwil3 = "Piwil3 significant"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff", Piwil3 = "black")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1, Piwil3 = 1)) +
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
                            str_c("all_biotype",
                                  "MA_plot", "DESeq2", "Mov10l1_KO_vs_WT",
                                  "Piwil3_KO_genes_highlight", "RefSeq",
                                  "png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)


### add gene names as labels
# filter data
plot_df_labels <-
  plot_df %>%
  dplyr::filter(regulation != "no") %>% 
  dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

# add labels
ma_plot_labeled <-
  ma_plot +
  geom_text(data = plot_df_labels,
            aes(x = mean, y = lfc, label = gene_name),
            check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
            colour = "black", fontface = "plain") +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("gray50", "red2", "#1a75ff", "black"))),
         alpha = F) +
  xlab("Mean expression") +
  ylab(str_c("log2 fold change: ", "Mov10l1 KO", " / ", "Mov10l1 WT", "\n")) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13, angle = 90)) +
  theme(legend.position = "bottom")

# save plot
ggsave(filename = file.path(outpath,
                            str_c("all_biotype",
                                  "MA_plot", "DESeq2", "Mov10l1_KO_vs_WT",
                                  "Piwil3_KO_genes_highlight", "RefSeq",
                                  "labeled",
                                  "png", sep = ".")),
       plot = ma_plot_labeled, width = 10, height = 10)
