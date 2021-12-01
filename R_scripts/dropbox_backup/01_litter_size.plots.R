### INFO: boxplot + jitter of litter size in two different genotype mice
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Zuzka")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
litter_df <- read_csv("Book3_ZL_littersize_LacZKO.csv")

######################################################## MAIN CODE
# clean data
litter_df_clean <- 
  litter_df %>% 
  dplyr::select(1:2) %>% 
  tidyr::gather(genotype, litter_size) %>% 
  dplyr::filter(!is.na(litter_size), 
                litter_size != 0)

# data statistics
litter_statistics <- 
  litter_df_clean %>% 
  group_by(genotype) %>% 
  summarize(median = median(litter_size), 
            sd = sd(litter_size), 
            length = length(litter_size)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

###### plots
# median + SD
set.seed(1000)
ggplot() +
  geom_jitter(data = litter_df_clean, aes(x = genotype, y = litter_size, fill = genotype), 
              shape = 21, colour = "black", size = 3, width = 0.2, height = 0.1, show.legend = F) +
  stat_summary(data = litter_df_clean, aes(x = genotype, y = litter_size), 
               fun.y = median, fun.ymin = median, fun.ymax = median, 
               geom = "crossbar", width = 0.5) +
  geom_errorbar(data = litter_statistics, aes(x = genotype, ymin = median - sd, ymax = median + sd), width = 0.3) +
  scale_fill_manual(values = c("red", "blue")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 15), breaks = 0:15) +
  xlab("genotype") +
  ylab("litter size") +
  ggsave(filename = "litter_size.median.sd.2.pdf", width = 10, height = 10)

# # boxplot
# ggplot(data = litter_df_clean, aes(x = genotype, y = litter_size)) +
#   geom_boxplot(show.legend = F, outlier.shape = NA, varwidth = T) +
#   geom_jitter(aes(colour = genotype), size = 3, width = 0.1, show.legend = F) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_y_continuous(breaks = 0:15) +
#   xlab("genotype") +
#   ylab("litter size") +
#   ggsave(filename = "litter_size.boxplot.pdf", width = 10, height = 10)


# # SEM only
# ggplot(data = litter_df_clean, aes(x = genotype, y = litter_size)) +
#   geom_jitter(aes(colour = genotype), size = 3, width = 0.1, show.legend = F) +
#   stat_summary(fun.y = median, 
#                fun.ymin = function(x) median(x) - (sd(x) / sqrt(length(x))), 
#                fun.ymax = function(x) median(x) + (sd(x) / sqrt(length(x))), 
#                geom = "errorbar", width = 0.5) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_y_continuous(breaks = 0:15) +
#   xlab("genotype") +
#   ylab("litter size") +
#   ggsave(filename = "litter_size.SEM.pdf", width = 10, height = 10)

# # median + SEM in box
# ggplot(data = litter_df_clean, aes(x = genotype, y = litter_size)) +
#   geom_jitter(aes(colour = genotype), size = 3, width = 0.1, show.legend = F) +
#   stat_summary(fun.y = median, 
#                fun.ymin = function(x) median(x) - (sd(x) / sqrt(length(x))), 
#                fun.ymax = function(x) median(x) + (sd(x) / sqrt(length(x))), 
#                geom = "crossbar", width = 0.5) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_y_continuous(breaks = 0:15) +
#   xlab("genotype") +
#   ylab("litter size") +
#   ggsave(filename = "litter_size.median.SEM.box.pdf", width = 10, height = 10)

# # median + SEM without box
# ggplot() +
#   geom_jitter(data = litter_df_clean, aes(x = genotype, y = litter_size, fill = genotype), 
#               shape = 21, colour = "black", size = 3, width = 0.3, height = 0, show.legend = F) +
#   stat_summary(data = litter_df_clean, aes(x = genotype, y = litter_size), 
#                fun.y = median, fun.ymin = median, fun.ymax = median, 
#                geom = "crossbar", width = 0.5) +
#   geom_errorbar(data = litter_statistics, aes(x = genotype, ymin = median - SEM, ymax = median + SEM), width = 0.3) +
#   scale_fill_manual(values = c("red", "blue")) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_y_continuous(limits = c(0, 15), breaks = 0:15) +
#   xlab("genotype") +
#   ylab("litter size") +
#   ggsave(filename = "litter_size.median.SEM.pdf", width = 10, height = 10)

# # median only
# ggplot() +
#   geom_jitter(data = litter_df_clean, aes(x = genotype, y = litter_size, fill = genotype), 
#               shape = 21, colour = "black", size = 3, width = 0.3, height = 0, show.legend = F) +
#   stat_summary(data = litter_df_clean, aes(x = genotype, y = litter_size), 
#                fun.y = median, fun.ymin = median, fun.ymax = median, 
#                geom = "crossbar", width = 0.5) +
#   scale_fill_manual(values = c("red", "blue")) + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_y_continuous(limits = c(0, 15), breaks = 0:15) +
#   xlab("genotype") +
#   ylab("litter size") +
#   ggsave(filename = "litter_size.median.pdf", width = 10, height = 10)

