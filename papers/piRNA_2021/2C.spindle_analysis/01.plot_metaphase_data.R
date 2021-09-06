### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/Svoboda/piRNA.Zuzka/Analysis/2020_paper/review/spindle_analysis")

######################################################## LIBRARIES
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(purrr)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# data path
data_path <- file.path(inpath, "MII_spindle_analysis_only.clean.xlsx")

######################################################## READ DATA
# read  
data_tb <- 
  openxlsx::read.xlsx(data_path, sheet = 1) %>% 
  as_tibble(.)

######################################################## MAIN CODE
# plot
purrr::map(names(data_tb)[1:3], function(name){
  
  # plot first two sheets with the same scale
  plot_tb <- 
    data_tb %>% 
    dplyr::select(animal, value = all_of(name)) %>% 
    dplyr::mutate(animal = str_extract(animal, "WT|KO"), 
                  animal = factor(animal, levels = c("WT", "KO")))
  
  # calculate mean and SD
  plot_tb_stats <- 
    plot_tb %>% 
    dplyr::group_by(animal) %>% 
    dplyr::summarise(N = n(), 
                     mean = mean(value), 
                     median = median(value),
                     sd = sd(value), 
                     stand_err = sd / sqrt(N)) %>% 
    dplyr::mutate(confidence_interval = stand_err * qt((0.95 / 2) + 0.5, N - 1)) 
  
  # plot
  bar_plot <- 
    ggplot() +
    geom_bar(data = plot_tb_stats, 
             aes(x = animal, y = mean, fill = animal),
             stat = "identity", 
             width = 0.85) +
    geom_errorbar(data = plot_tb_stats, 
                  aes(x = animal, y = mean, ymin = mean - sd, ymax = mean + sd), 
                  color = "gray50",
                  size = 3, 
                  width = 0.45) +
    geom_errorbar(data = plot_tb_stats, 
                  aes(x = animal, y = mean, ymin = mean, ymax = mean), 
                  color = "gray50",
                  size = 3,
                  width = 0.85) +
    geom_jitter(data = plot_tb,
                mapping = aes(x = animal, y = value, fill = animal),
                color = "black", 
                size = 16, 
                shape = 15, 
                position = position_jitterdodge(dodge.width = 1, jitter.width = 0.6)) +
    scale_fill_manual(values = c("WT" = "gray75", "KO" = "gray25")) +
    scale_color_manual(values = c("WT" = "black", "KO" = "black")) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), 
          axis.text.y = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none")
  
  # save
  ggsave(plot = bar_plot, 
         filename = file.path(outpath, str_c(str_c(name, "_barplot"), "mean", "standard_deviation", "wmf", sep = ".")),
         width = 8,
         height = 20)
  
  # return
  return(name)
  
})


