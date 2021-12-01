### INFO: 
### DATE: Fri Aug 07 19:24:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.MesAur1/LTRs/LTR_exaptations.GenRes")

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

# exaptated genes path
exaptated_genes_path <- file.path(inpath, "supp_gr.216150.116_Supplemental_Table_S5.xlsx")

# diff. exp. results path
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression"
results_path <- list.files(results_path, ".*\\.all_results\\.xlsx", recursive = T, full.names = T)

######################################################## READ DATA
# read exaptated genes
exaptated_genes <- 
  openxlsx::read.xlsx(exaptated_genes_path, startRow = 2) %>% 
  as_tibble(.)

# read diff. exp. results
results_tb <- 
  openxlsx::read.xlsx(results_path, sheet = "Mov10l_KO_vs_Mov10l_WT") %>% 
  as_tibble(.)

######################################################## MAIN CODE
# clean exaptated genes
exaptated_genes_clean <- 
  exaptated_genes %>% 
  dplyr::select(gene_name = host.gene, exaptation.type, LTR.group)

## MA plot
# get axis limits
results_limits <-
  results_tb %>%
  dplyr::summarise(x_limit = baseMean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = log2FoldChange %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# data for plot
plot_df <-
  results_tb %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name) %>%
  dplyr::left_join(., exaptated_genes_clean, by = "gene_name") %>%
  dplyr::mutate(exaptation.type = replace(exaptation.type, is.na(exaptation.type), "no exaptation"), 
                exaptation.type = factor(exaptation.type, levels = c("5' exon", "3' exon", "internal exon", "tagging", "no exaptation"))) %>%
  dplyr::group_by(gene_name) %>% 
  dplyr::slice_sample(n = 1) %>% 
  dplyr::arrange(desc(exaptation.type))
  
# plot
ma_plot <-
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = exaptation.type, alpha = exaptation.type), size = 2.5, shape = 20) +
  scale_x_log10(limits = c(0.01, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                     breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
  scale_colour_manual(values = c(`no exaptation` = "gray50", 
                                 `5' exon` = "forestgreen", 
                                 `3' exon` = "gray50", 
                                 `internal exon` = "gray50",
                                 tagging = "gray50")) +
  scale_alpha_manual(values = c(`no exaptation` = 0.5, 
                                `5' exon` = 1, 
                                `3' exon` = 1, 
                                `internal exon` = 1,
                                tagging = 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

### add gene names as labels
# # filter data
# plot_df_labels <-
#   plot_df %>%
#   dplyr::filter(exaptation.type != "no exaptation") %>%
#   dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

# add labels
ma_plot_labeled <-
  ma_plot +
  # geom_text(data = plot_df_labels,
  #           aes(x = mean, y = lfc, label = gene_name),
  #           check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
  #           colour = "black", fontface = "plain") +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = rev(c("gray50", "coral", "red2", "royalblue", "forestgreen")))),
         alpha = F) +
  xlab("Mean expression") +
  ylab(str_c("log2 fold change: ", "Mov10l1 KO", " / ", "Mov10l1 WT", "\n") %>% str_replace_all(., "_", " ")) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13, angle = 90)) +
  theme(legend.position = "bottom")

# save plot
ggsave(filename = file.path(outpath, str_c("Mov10l_KO_vs_Mov10l_WT.oocyte.GenRes.exapatated_genes.type.png", sep = ".")),
       plot = ma_plot_labeled, width = 10, height = 10)

