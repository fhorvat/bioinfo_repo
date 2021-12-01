#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: counts reads in categories (rRNA, repeat, exon, other) in chunks
### DATE: Mon Mar 04 14:51:59 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/datasets/mouse_testis.Papd7.small_RNAseq.Apr_2021/Analysis/read_class")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(data.table)

library(cowplot)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# list read class files
read_class_path <- list.files(inpath, "\\.read_class\\.widths\\.csv", full.names = T, recursive = T)

######################################################## READ DATA
# read class files
read_class_list <- 
  purrr::map(read_class_path, readr::read_csv) %>% 
  dplyr::bind_rows(.) 

######################################################## MAIN CODE
# set feature classes
class_hier <- 
  read_class_list %>% 
  dplyr::select(read_class = gene_biotype, read_class_precise = read_group) %>% 
  dplyr::distinct(.) %>% 
  dplyr::mutate(read_class = ifelse(is.na(read_class), read_class_precise, read_class),
                read_class = str_replace_all(read_class, c("ENSEMBL" = "Ens.", "RepeatMasker" = "Rmsk.")),
                read_class = factor(read_class, levels = rev(c("IAP", "LTR", "SINE", "LINE", 
                                                               "DNA/simple repeat", "non-coding RNA Rmsk.",
                                                               "other annotated Rmsk.", 
                                                               "non-coding RNA Ens.", 
                                                               "other annotated Ens.", "protein coding Ens.", 
                                                               "rRNA", "tRNA", "miRNA",
                                                               "not annotated"))), 
                read_class_precise = factor(read_class_precise, unique(read_class_precise))) 

# filter table
read_class_tb <- 
  read_class_list %>% 
  dplyr::select(read_width, count, read_class = gene_biotype, read_class_precise = read_group, sample_id) %>% 
  dplyr::mutate(read_class = str_replace_all(read_class, c("ENSEMBL" = "Ens.", "RepeatMasker" = "Rmsk."))) %>% 
  dplyr::mutate(read_class_precise = factor(read_class_precise, unique(read_class_precise))) %>% 
  dplyr::filter((read_width >= 25) & (read_width <= 32)) %>% 
  dplyr::group_by(sample_id, read_class_precise) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(count_fract = (count / total_count)) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Papd7_WT|Papd7_HET|Papd7_KO|Papd7_ex4_KO") %>% 
                  str_replace(., "Papd7_", "Papd7 "), 
                age = str_extract(sample_id, "7dpp|14dpp")) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Papd7 WT", "Papd7 HET", "Papd7 KO", "Papd7 ex4_KO")), 
                age = factor(age, levels = c("7dpp", "14dpp")))

# mean values for age/genotype
read_class_mean <- 
  read_class_tb %>% 
  dplyr::group_by(age, genotype, read_class_precise) %>% 
  dplyr::summarise(count_fract = mean(count_fract)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(age)



### plot for individual samples
# get data for plot
plot_tb <- 
  read_class_tb %>% 
  dplyr::left_join(., class_hier, by = "read_class_precise") %>% 
  dplyr::group_by(sample_id, genotype, age, read_class) %>% 
  dplyr::summarise(count_fract = sum(count_fract)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::arrange(age, genotype) %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "^s_testis_|_r[0-9]+\\.SE.*$") %>% str_replace_all(., "_", " "), 
                sample_id = factor(sample_id, levels = unique(sample_id))) %>% 
  dplyr::arrange(sample_id) 

# plot
barplot_individual <- 
  ggplot() + 
  geom_bar(data = plot_tb, 
           mapping = aes(x = sample_id, y = count_fract, fill = read_class), 
           colour = "black", width = 0.9, stat = "identity") + 
  scale_fill_manual(values = c("white",
                               RColorBrewer::brewer.pal(3, "Reds"), 
                               RColorBrewer::brewer.pal(8, "Purples")[3:8], 
                               RColorBrewer::brewer.pal(4, "Blues"))) + 
  guides(fill = guide_legend(reverse = FALSE)) +
  labs(subtitle = "Anontation of mapped small RNAs (25-32 nt) - individual samples") +
  xlab("Sample name") +
  ylab("Ratio") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.subtitle = element_text(size = 18), 
        legend.position = "right", 
        # legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 



### plot for individual samples precise TE classes
# get data for plot
plot_tb <- 
  read_class_list %>% 
  dplyr::select(read_width, count, read_class = gene_biotype, read_class_precise = read_group, sample_id) %>% 
  dplyr::filter(read_class %in% c("SINE", "LINE", "LTR")) %>% 
  dplyr::mutate(read_class_precise = factor(read_class_precise, unique(read_class_precise))) %>% 
  dplyr::filter((read_width >= 25) & (read_width <= 32)) %>% 
  dplyr::group_by(sample_id, read_class_precise) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(count_fract = (count / total_count)) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Papd7_WT|Papd7_HET|Papd7_KO|Papd7_ex4_KO") %>% 
                  str_replace(., "Papd7_", "Papd7 "), 
                age = str_extract(sample_id, "7dpp|14dpp")) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Papd7 WT", "Papd7 HET", "Papd7 KO", "Papd7 ex4_KO")), 
                age = factor(age, levels = c("7dpp", "14dpp"))) %>% 
  dplyr::arrange(age, genotype) %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "^s_testis_|_r[0-9]+\\.SE.*$") %>% str_replace_all(., "_", " "), 
                sample_id = factor(sample_id, levels = unique(sample_id))) %>% 
  dplyr::arrange(sample_id) 

# plot
barplot_individual_precise <- 
  ggplot() + 
  geom_bar(data = plot_tb, 
           mapping = aes(x = sample_id, y = count_fract, fill = read_class_precise), 
           colour = "black", width = 0.9, stat = "identity") + 
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Blues")[4:8], 
                               RColorBrewer::brewer.pal(8, "PuRd")[3:8], 
                               RColorBrewer::brewer.pal(8, "Greys")[2:8])) + 
  guides(fill = guide_legend(reverse = FALSE)) +
  labs(subtitle = "Anontation of TE annotated small RNAs (25-32 nt) - individual samples") +
  xlab("Sample name") +
  ylab("Ratio") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.subtitle = element_text(size = 18), 
        legend.position = "right", 
        # legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) 



### plot for each age as facet
# get data for plot
plot_tb <- 
  read_class_mean %>% 
  dplyr::left_join(., class_hier, by = "read_class_precise") %>% 
  dplyr::group_by(genotype, age, read_class) %>% 
  dplyr::summarise(count_fract = sum(count_fract)) %>% 
  dplyr::ungroup(.)

# plot
barplot_mean <- 
  ggplot() + 
  geom_bar(data = plot_tb, 
           mapping = aes(x = genotype, y = count_fract, fill = read_class), 
           colour = "black", width = 0.9, stat = "identity") + 
  facet_grid(cols = vars(age), scales = "free", drop = T) + 
  scale_fill_manual(values = c("white",
                               RColorBrewer::brewer.pal(3, "Reds"), 
                               RColorBrewer::brewer.pal(8, "Purples")[3:8], 
                               RColorBrewer::brewer.pal(4, "Blues"))) + 
  guides(fill = guide_legend(reverse = FALSE)) +
  labs(subtitle = "Anontation of mapped small RNAs (25-32 nt) - mean per genotype") +
  xlab("Genotype") +
  ylab("Ratio") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.subtitle = element_text(size = 18), 
        legend.position = "right", 
        # legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  theme(strip.text.x = element_text(size = 16))


### plot  precise for LTRs, LINEs and SINEs subfamilies
# get data for plot
plot_tb <- 
  read_class_list %>% 
  dplyr::select(read_width, count, read_class = gene_biotype, read_class_precise = read_group, sample_id) %>% 
  dplyr::filter(read_class %in% c("SINE", "LINE", "LTR")) %>% 
  dplyr::mutate(read_class_precise = factor(read_class_precise, unique(read_class_precise))) %>% 
  dplyr::filter((read_width >= 25) & (read_width <= 32)) %>% 
  dplyr::group_by(sample_id, read_class_precise) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(total_count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(count_fract = (count / total_count)) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Papd7_WT|Papd7_HET|Papd7_KO|Papd7_ex4_KO") %>% 
                  str_replace(., "Papd7_", "Papd7 "), 
                age = str_extract(sample_id, "7dpp|14dpp")) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Papd7 WT", "Papd7 HET", "Papd7 KO", "Papd7 ex4_KO")), 
                age = factor(age, levels = c("7dpp", "14dpp"))) %>% 
  dplyr::group_by(age, genotype, read_class_precise) %>% 
  dplyr::summarise(count_fract = mean(count_fract)) %>% 
  dplyr::ungroup(.)

# plot
barplot_mean_precise <- 
  ggplot() + 
  geom_bar(data = plot_tb, 
           mapping = aes(x = genotype, y = count_fract, fill = read_class_precise), 
           colour = "black", width = 0.9, stat = "identity") + 
  facet_grid(cols = vars(age), scales = "free", drop = T) + 
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Blues")[4:8], 
                               RColorBrewer::brewer.pal(8, "PuRd")[3:8], 
                               RColorBrewer::brewer.pal(8, "Greys")[2:8])) + 
  guides(fill = guide_legend(reverse = FALSE)) +
  labs(subtitle = "Anontation of TE annotated small RNAs (25-32 nt) - mean per genotype") +
  xlab("Genotype") +
  ylab("Ratio") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.subtitle = element_text(size = 18), 
        legend.position = "right", 
        # legend.key.size = unit(1, "cm"), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  theme(strip.text.x = element_text(size = 16))


### grid and save
# grid 
grid <- plot_grid(barplot_individual, barplot_individual_precise, 
                  barplot_mean, barplot_mean_precise,
                  ncol = 2,
                  nrow = 2,
                  align = "v", 
                  axis = "b") 

# # save as pdf
# ggsave(filename = file.path(outpath, "read_class_distribution.pdf"),
#        plot = grid,
#        width = 20,
#        height = 20)

# save as object
saveRDS(grid, file.path(outpath, "grid_2.RDS"))


