### INFO: 
### DATE: Tue Aug 03 14:21:57 2021
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/functional_oocyte_miRNA/maternal_miRNA_expression/expression/cow.bosTau9")

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

# expression paths
fpm_path <- list.files(inpath, ".*\\.miRNA\\.FPM\\.csv", full.names = T)
fpm_mean_path <- list.files(inpath, ".*\\.miRNA\\.FPM_mean\\.csv", full.names = T)

######################################################## READ DATA
# read expression
fpm_tb <- readr::read_csv(fpm_path)
fpm_mean <- readr::read_csv(fpm_mean_path)

######################################################## MAIN CODE
# get the top 100 most abundant miRNAs
mirna_top <- 
  fpm_mean %>% 
  dplyr::top_n(n = 100, wt = GV) %>% 
  dplyr::arrange(-GV) %$% 
  gene_id
  
# filter table
fpm_tb_top <- 
  fpm_tb %>% 
  dplyr::filter(gene_id %in% mirna_top) %>% 
  dplyr::arrange(match(gene_id, mirna_top))

# get long table
plot_tb <- 
  fpm_tb_top %>% 
  dplyr::select(-c(gene_biotype, gene_description)) %>% 
  tidyr::pivot_longer(cols = -c(gene_id, gene_name, coordinates, strand), names_to = "sample_id", values_to = "FPM") %>% 
  dplyr::mutate(stage = str_extract(sample_id, "GV|MI(?=_)|MII(?=_)"), 
                stage = factor(stage, levels = c("GV", "MI", "MII"))) %>% 
  dplyr::arrange(stage) %>% 
  dplyr::mutate(sample_id = factor(sample_id, levels = unique(sample_id)))
  
# plot as heatmap ggplot2
heat_plot <- 
  ggplot(plot_tb, aes(x = sample_id, y = gene_id)) + 
  geom_tile(aes(fill = FPM), colour = "grey45") + 
  # coord_equal() + 
  scale_fill_gradient2(low = "gray95", mid = "cyan", high = "darkgreen") +
  theme(axis.text.x = element_text(size = 12, angle = 90, face = "bold", colour = "grey25", vjust = 0.5, hjust = 0), 
        axis.text.y = element_text(size = 12, face = "bold", colour = "grey25"), 
        legend.title = element_text(size = 10, face = "bold"), 
        # legend.position = "rigth", 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = NA), 
        axis.ticks = element_blank()) + 
  labs(x = "", 
       y = "", 
       fill = "") + 
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = levels(plot_tb$gene_id), position = "right")

# save
ggsave(plot = heat_plot, filename = file.path(outpath, str_c("test", ".20210225.png")), width = 15, height = 6)
