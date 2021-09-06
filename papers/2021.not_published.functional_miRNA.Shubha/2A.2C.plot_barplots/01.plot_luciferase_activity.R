### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("C:/Users/Petr/Dropbox/Bioinfo/Svoboda/miRNA.Shubha/functional_oocyte_miRNA/luciferase_acivity_plots")

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
data_path <- file.path(inpath, "210608 For Filip (luciferase).xlsx")

######################################################## READ DATA
# get sheet names
sheet_names <- openxlsx::getSheetNames(data_path)
  
# read sheets
data_tb <- purrr::map(sheet_names, function(sheet){

  # read 
  data_tb <- 
    openxlsx::read.xlsx(data_path, sheet) %>% 
    as_tibble(.) %>% 
    dplyr::select(1:3) %>% 
    dplyr::mutate(experiment = str_remove(sheet, "\\(.*"))
  
  # return
  return(data_tb)
  
}) %>%
  set_names(str_remove(sheet_names, "\\(.*"))

######################################################## MAIN CODE
### plot first two sheets with the same scale
# get the y-scale
y_scale <- 
  data_tb[1:2] %>% 
  dplyr::bind_rows(.) %$% 
  Relative.nanoluc.activity %>% 
  max(.)

# set the y-scale
y_scale <- 1.5

# plot
purrr::map(names(data_tb)[1:2], function(name){
  
  # plot first two sheets with the same scale
  plot_tb <- 
    data_tb[[name]] %>% 
    dplyr::mutate(Reporter = factor(Reporter, levels = c("1x-perfect", "4x-bulged", "4x-mutant")))
  
  # calculate mean and SD
  plot_tb_stats <- 
    plot_tb %>% 
    dplyr::group_by(Reporter) %>% 
    dplyr::summarise(N = n(), 
                     mean = mean(Relative.nanoluc.activity), 
                     median = median(Relative.nanoluc.activity),
                     sd = sd(Relative.nanoluc.activity), 
                     stand_err = sd / sqrt(N)) %>% 
    dplyr::mutate(confidence_interval = stand_err * qt((0.95 / 2) + 0.5, N - 1)) 
  
  # plot
  bar_plot <- 
    ggplot() +
    geom_bar(data = plot_tb_stats, 
             aes(x = Reporter, y = mean),
             fill = "#d9d9d9",
             stat = "identity", 
             width = 0.95) +
    geom_errorbar(data = plot_tb_stats, 
                  aes(x = Reporter, y = mean, ymin = mean - sd, ymax = mean + sd), 
                  color = "#7f7f7f",
                  size = 3, 
                  width = 0.55) +
    geom_errorbar(data = plot_tb_stats, 
                  aes(x = Reporter, y = mean, ymin = mean, ymax = mean), 
                  color = "#404040",
                  size = 3,
                  width = 0.95) +
    geom_jitter(data = plot_tb,
                mapping = aes(x = Reporter, y = Relative.nanoluc.activity, color = Reporter),
                size = 8, 
                shape = 15, 
                position = position_jitterdodge(dodge.width = 1, jitter.width = 0.6)) +
    scale_color_manual(values = c("black", "black", "black")) +
    scale_y_continuous(limits = c(0, y_scale), 
                       breaks = seq(0, y_scale, 0.1)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), 
          axis.text.y = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.position = "none")  +
  ggsave(filename = file.path(outpath, str_c("luciferase_barplot", "mean", "standard_deviation", name, "wmf", sep = ".")),
         width = 8,
         height = 20,
         device = "wmf")

  # return
  return(name)
  
})


### plot the last sheet separately
# set name
name <- names(data_tb)[3]

# plot first two sheets with the same scale
plot_tb <- 
  data_tb[[name]] %>% 
  dplyr::rename(blast_perc = `%.blastocysts`) %>% 
  dplyr::mutate(Injected = stringr::str_squish(Injected),
                Injected = factor(Injected, levels = c("water", "Let-7a antagomir", "miR-205 antagomir")))

# calculate mean and SD
plot_tb_stats <- 
  plot_tb %>% 
  dplyr::group_by(Injected) %>% 
  dplyr::summarise(N = n(), 
                   mean = mean(blast_perc), 
                   median = median(blast_perc),
                   sd = sd(blast_perc), 
                   stand_err = sd / sqrt(N)) %>% 
  dplyr::mutate(confidence_interval = stand_err * qt((0.95 / 2) + 0.5, N - 1)) 

# set the y-scale
y_scale <- 20

# plot
bar_plot <- 
  ggplot() +
  geom_bar(data = plot_tb_stats, 
           aes(x = Injected, y = mean),
           fill = "#d9d9d9",
           stat = "identity", 
           width = 0.95) +
  geom_errorbar(data = plot_tb_stats, 
                aes(x = Injected, y = mean, ymin = mean - sd, ymax = mean + sd), 
                color = "#7f7f7f",
                size = 3, 
                width = 0.55) +
  geom_errorbar(data = plot_tb_stats, 
                aes(x = Injected, y = mean, ymin = mean, ymax = mean), 
                color = "#404040",
                size = 3,
                width = 0.95) +
  geom_jitter(data = plot_tb,
              mapping = aes(x = Injected, y = blast_perc, color = Injected),
              size = 8, 
              shape = 15, 
              position = position_jitterdodge(dodge.width = 1, jitter.width = 0.6)) +
  scale_color_manual(values = c("black", "black", "black")) +
  scale_y_continuous(limits = c(0, y_scale), 
                     breaks = seq(0, y_scale, 2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20, vjust = 0.5), 
        axis.text.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none") +
  ggsave(filename = file.path(outpath, str_c("luciferase_barplot", "mean", "standard_deviation", name, "wmf", sep = ".")), 
         width = 8, 
         height = 20,
         device = "wmf")
