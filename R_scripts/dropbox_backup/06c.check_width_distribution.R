### INFO: 
### DATE: Fri Dec 11 17:05:38 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/testis.RNA_seq/distance_to_retrotransposons")

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

library(GenomicRanges)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets"

# experiment
experiment_name <- "hamster_testis_Mov10l.8.5dpp.run_2.RNAseq"
experiment_name <- "hamster_testis_Mov10l.0.5dpp.RNAseq"

# results path
expression_path <- file.path(base_path, experiment_name,
                             "Analysis/expression.added_PIWIL3.stranded", 
                             "results.ensembl.99.MesAur1.0.20200415.UCSCseqnames.Piwil3_Tex101.from_RefSeq")
results_path <- file.path(expression_path, "protein_coding.diffExp.DESeq2.genotype_age.significant_results.xlsx")
results_path_all <- file.path(expression_path, "protein_coding.diffExp.DESeq2.genotype_age.all_results.xlsx")

######################################################## READ DATA
# read differently expressed genes
results_tb <- openxlsx::read.xlsx(results_path) %>% as_tibble(.)

# read all genes
results_tb_all <- openxlsx::read.xlsx(results_path_all) %>% as_tibble(.)

######################################################## MAIN CODE
### tidy data
# get GRanges of all genes
results_gr_all <- 
  results_tb_all %>% 
  dplyr::select(gene_id, log2FoldChange, baseMean, coordinates) %>%
  tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.)


### plot relationship between logFC and gene length
# join all
results_tb_all <- 
  results_gr_all %>% 
  as_tibble(.) %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down"), 
                regulation = replace(regulation, !(gene_id %in% results_tb$gene_id), "no"), 
                regulation = factor(regulation, levels = c("no", "up", "down"))) %>% 
  dplyr::arrange(regulation)


axis_limits <- 4

# scatter plot
scatter_plot <-
  ggplot(results_tb_all, aes(x = width, y = log2FoldChange, color = regulation, alpha = regulation)) +
  geom_point(shape = 16, size = 3) +
  scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                      values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  # scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  xlab(str_c("gene length")) +
  ylab(str_c("log2FC DESeq2 Mov10l1 KO vs. Mov10l1 WT")) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("scatter_plot", "log2FC_vs_gene_length", experiment_name, "png", sep = ".")),
       plot = scatter_plot,
       width = 10, height = 10)



# # get GRanges of upregulated genes
# results_gr <- 
#   results_tb %>% 
#   dplyr::filter(log2FoldChange > 0) %>%
#   dplyr::select(gene_id, log2FoldChange, baseMean, coordinates) %>%
#   tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
#   GRanges(.)
# 
# ### get set of random genes
# # repeat 1000 times
# set.seed(1234)
# width_random <- purrr::map(1:1000, function(n){
#   
#   # sample random genes
#   random_gr <- results_gr_all[sample(1:length(results_gr_all), length(results_gr), replace = F)]
#   
#   # get mean distance
#   return(mean(width(random_gr)))
#   
# })
# 
# # find quantiles
# width_tb <-
#   tibble(average_width = c(mean(width(results_gr)),
#                            unlist(width_random)),
#          source = c(experiment_name,
#                     rep("random", length(width_random)))) %>%
#   dplyr::mutate(promile = ntile(average_width, 1000)) %>%
#   dplyr::mutate(percentile = ntile(average_width, 100)) %>%
#   dplyr::arrange(-average_width) %>%
#   dplyr::filter(source == experiment_name)




