### INFO: plots Kimura's divergence for repeats (clone of createRepeatLandscape.pl)
### DATE: Thu Jan 20 07:53:06 2022
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
div_path <- args$div_path
genome_2bit_path <- args$genome_2bit_path
classes_path <- args$classes_path
colors_path <- args$colors_path
plot_name <- args$plot_name

######################################################## READ DATA
# get genome length (without Ns)
genome_lengths <- system(command = str_c("twoBitInfo -noNs ", genome_2bit_path, " stdout | awk '{sum+=$2}END{print sum}'"), 
                         intern = T) %>% 
  as.numeric(.)

# read divergence file
div_tb <- data.table::fread(div_path, skip = "Coverage for each repeat")

# read classes table
classes_tb <- readr::read_delim(classes_path, delim = "\t", col_names = c("repClass.Family", "repeat_class"))

# read colors table
colors_tb <- readr::read_delim(colors_path, delim = "\t", col_names = c("repeat_class", "rep_color"))

######################################################## MAIN CODE
# clean table, join with categories, calculate stuff
div_tb_clean <-
  div_tb %>% 
  tidyr::pivot_longer(cols = -Div, names_to = "repClass.Family", values_to = "length") %>% 
  dplyr::left_join(., classes_tb, by = "repClass.Family") %>% 
  dplyr::mutate(repeat_class = ifelse(is.na(repeat_class), repClass.Family, repeat_class)) %>% 
  dplyr::left_join(., colors_tb, by = "repeat_class") %>% 
  dplyr::filter(!is.na(rep_color)) %>% 
  dplyr::mutate(length_perc = 100*(length / genome_lengths)) %>% 
  dplyr::rename(kimura_div = Div)

# create table for plot
plot_tb <- 
  div_tb_clean %>% 
  dplyr::filter(!(repeat_class %in% c("Simple_repeat", "Satellite", "Structural_RNA"))) %>% 
  dplyr::group_by(repeat_class, kimura_div) %>% 
  dplyr::summarise(length_perc = sum(length_perc)) %>%  
  dplyr::ungroup(.) %>% 
  dplyr::mutate(repeat_class = factor(repeat_class, levels = rev(unique(colors_tb$repeat_class[colors_tb$repeat_class %in% repeat_class])))) %>% 
  dplyr::arrange(repeat_class)

# get y-limit
y_limit <- 
  plot_tb %>% 
  dplyr::group_by(kimura_div) %>% 
  dplyr::summarise(length_perc = sum(length_perc)) %$% 
  length_perc %>% 
  .[which.max(.)]

# round y-limit 
y_limit <- round(y_limit + 10 * 10^-2, 1)

# create color scale
color_scale <- colors_tb$rep_color
names(color_scale) <- colors_tb$repeat_class
color_scale <- color_scale[names(color_scale) %in% plot_tb$repeat_class]
color_scale <- c(rev(color_scale[str_detect(names(color_scale), "^SINE")]), 
                 rev(color_scale[str_detect(names(color_scale), "^LINE")]), 
                 rev(color_scale[str_detect(names(color_scale), "^LTR")]), 
                 color_scale[names(color_scale) == "RC/Helitron"], 
                 rev(color_scale[str_detect(names(color_scale), "^DNA")]),
                 color_scale[names(color_scale) == "Unknown"])

# plot for each sample
kimura_barplot <- 
  ggplot() + 
  geom_bar(data = plot_tb, 
           mapping = aes(x = kimura_div, y = length_perc, fill = repeat_class), 
           width = 0.8, stat = "identity", color = "white") + 
  scale_x_continuous(limits = c(0, 51), breaks = seq(0, 50, 5)) +
  scale_y_continuous(limits = c(0, y_limit), breaks = seq(0, y_limit, 0.2)) + 
  coord_cartesian(xlim = c(2.1, 48.4), 
                  ylim = c(0.085, y_limit)) + 
  scale_fill_manual(values = color_scale) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, face = "italic", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "italic", margin = margin(t = 0, r = 20, b = 0, l = 10)),
        axis.text.x = element_text(size = 15, vjust = 0.5, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 15, face = "bold", margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        legend.margin = margin(t = 0, r = 20, b = 0, l = 20), 
        axis.ticks = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.border = element_blank(),
        legend.position = "right", 
        legend.title = element_blank()) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = "Kimura substitution level (CpG adjusted)",
       y = "percent of genome", 
       title = "Interspersed Repeat Landscape")

# save
ggsave(plot = kimura_barplot, 
       filename = file.path(outpath, str_c(plot_name, ".png")), 
       width = 18, height = 10)
