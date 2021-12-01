### INFO: boxplot + jitter of litter size in two different genotype mice
### DATE: 06. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/Zuzka/2018")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggbeeswarm)
library(ggplus)
library(gridExtra)
library(ggpubr)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## READ DATA
litter_df <- read_csv("Litter size data PAPD7.csv")

######################################################## MAIN CODE
# clean data
litter_df %<>% 
  dplyr::mutate_at(.vars = vars(dplyr::starts_with("litter")), .funs = funs(replace(., is.na(.), 0))) %>% 
  tidyr::gather(key = litter, value = litter_size, 3:ncol(.)) %>% 
  dplyr::select(graph, genotype, litter_size)

litter_df %>% 
  filter(graph == "graph1") %$% 
  litter_size %>% 
  table

# split by graphs
litter_df_list <- split(litter_df, f = litter_df$graph)

# calculate statistics, plot
invisible(lapply(litter_df_list, function(litter_df_clean){
  
  # set factors, remove 0 size litters
  litter_df_clean %<>%
    mutate(genotype = as.factor(genotype)) %>% 
    dplyr::filter(litter_size != 0)
  
  # set title
  title <- unique(litter_df_clean$graph)
  
  # data statistics
  litter_statistics <- 
    litter_df_clean %>% 
    group_by(genotype) %>% 
    summarize(median = median(litter_size), 
              sd = sd(litter_size), 
              length = length(litter_size)) %>% 
    dplyr::mutate(SEM = sd / sqrt(length))
  
  # litter comparison
  if(length(unique(litter_df_clean$genotype)) > 1){
    p_label <- wilcox.test(litter_df_clean %>% dplyr::filter(genotype == unique(litter_df_clean$genotype)[1]) %$% litter_size, 
                           litter_df_clean %>% dplyr::filter(genotype == unique(litter_df_clean$genotype)[2]) %$% litter_size)
    if(p_label$p.value <= 0.05){
      p_label <- str_c("Wilcoxon, p = ", round(p_label$p.value, digits = 6))
    }else{
      p_label <- "not significant"
    }
  }else{
    p_label <- "not significant"
  }
  
  # median + SD dotplot
  ggplot() +
    geom_dotplot(data = litter_df_clean, aes(x = genotype, y = litter_size, fill = genotype), 
                 binaxis = "y", binwidth = 1, dotsize	= 0.25, stackdir = "center", colour = "black", show.legend = F) +
    stat_summary(data = litter_df_clean, aes(x = genotype, y = litter_size), 
                 fun.y = median, fun.ymin = median, fun.ymax = median, 
                 geom = "crossbar", width = 0.3) +
    geom_errorbar(data = litter_statistics, aes(x = genotype, ymin = median - sd, ymax = median + sd), width = 0.2) +
    annotate("text", label = p_label, x = 0.6, y = 16, size = 3, colour = "black") + 
    scale_fill_manual(values = c("red", "blue")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(0, 16), breaks = 0:15) +
    ggtitle(title) +
    xlab("genotype") +
    ylab("litter size") +
    ggsave(filename = str_c(title, ".litter_size.dotplot.median.sd.pdf"), width = 10, height = 10)
  
}))
