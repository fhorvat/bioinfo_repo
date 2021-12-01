rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/final_figures/figure_6B_S5B_S5C_S5D")

######################################################## LIBRARIES
library(dplyr)
library(tidyr)
library(magrittr)
library(readr)
library(stringr)
library(ggplot2)

######################################################## FUNCTIONS
plotCummulativeFPKMBin <- function(coverage_df, tile_width, experiment_name, original_sample_order, plot_sample_order){

  # set color palette 
  color_palette <- c("black", "grey60", "orange", "grey30", "grey50")
  color_palette <- color_palette[1:length(original_sample_order)]
    
  # wide format
  coverage_wide <- 
    coverage_df %>% 
    dplyr::select(stage, pos, count, element_width, element_name) %>% 
    mutate(count = replace(count, is.na(count), 0)) %>% 
    tidyr::spread(key = stage, value = count) %>% 
    data.table::setnames(x = ., old = 4:ncol(.), new = str_c("bam_", colnames(.[4:ncol(.)])))

  # bin
  coverage_wide_binned <- 
    rbind(coverage_wide %>% #upstream
            filter(pos >= -150000 & pos < 0) %>%
            group_by(indx = gl(ceiling(n()/tile_width), tile_width, n())) %>%
            dplyr::select(4:ncol(.)) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = paste(seq(150000, tile_width, -tile_width), "upstream", sep = "_")),
          coverage_wide %>% #element
            filter(pos >= 0 & pos <= coverage_wide$element_width[1]) %>%
            dplyr::select(4:ncol(.)) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = "element", 
                   indx = 1) %>%
            dplyr::select(c(ncol(.), 1:(ncol(.) - 1))), 
          coverage_wide %>% #downstream
            filter(pos > coverage_wide$element_width[1] & pos <= (150000 + coverage_wide$element_width[1])) %>%
            group_by(indx = gl(ceiling(n()/tile_width), tile_width, n())) %>%
            dplyr::select(4:ncol(.)) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_"))) %>%
    dplyr::select(-1)
  
  # FPM to FPKM (for bins)
  coverage_wide_binned[, 1:(ncol(coverage_wide_binned) - 1)] <- coverage_wide_binned[, 1:(ncol(coverage_wide_binned) - 1)] / (tile_width / 1000)
  
  # output table
  coverage_wide_binned %>%
    mutate(pos = c(seq(-150000, -tile_width, tile_width), 0, seq(coverage_wide$element_width[1], ((150000 + coverage_wide$element_width[1]) - tile_width), tile_width))) %T>%
    write.csv(file = str_c(experiment_name, "_", coverage_wide$element_name[1], "_cummulativeFPKMBins_", tile_width / 1000, "kb.csv"), row.names = F)
  
  # get max y-lim
  max_ylim <- 
    coverage_wide_binned %>% 
    dplyr::filter(position != "element") %>% 
    dplyr::select(1:(ncol(.) - 1)) %>% 
    max(.)
  
  # long format
  coverage_long_binned <- 
    coverage_wide_binned %>% 
    tidyr::gather(key = stage, value = fpkm, -position) %>% 
    mutate(element = ifelse((position == "element"), "in_element",  "out_element"), 
           width = ifelse((position == "element"), coverage_wide$element_width[1], tile_width), 
           pos = rep(c(seq(-150000, -tile_width, tile_width), 0, 
                       seq(coverage_wide$element_width[1], ((150000 + coverage_wide$element_width[1]) - tile_width), tile_width)), 
                     length(original_sample_order)), 
           stage = str_replace_all(stage, "bam_", "")) %>%
    dplyr::select(pos, fpkm, element, stage, position, width) %>%
    mutate(stage = factor(stage, levels = plot_sample_order)) %>%
    arrange(stage) %>%
    mutate(stage = factor(stage, levels = original_sample_order))
  
  # plot
  ggplot(coverage_long_binned, aes(x = pos, y = fpkm)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + width, ymax = fpkm, fill = stage)) +
    scale_fill_manual(values = color_palette, 
                      breaks = original_sample_order) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + coverage_wide$element_width[1]), 
                       breaks = c(seq(-150000, -10000, 10000), 0, seq(coverage_wide$element_width[1], (140000 + coverage_wide$element_width[1]), 10000)), 
                       labels = c(seq(150, 10, -10), "element", seq(0, 140, 10)), 
                       name = "bin") + 
    scale_y_continuous(name = "Cummulative FPKM") +  
    coord_cartesian(ylim = c(0, max_ylim + 50)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = str_c(experiment_name, "_", coverage_wide$element_name[1], "_cummulativeFPKMBins_", tile_width / 1000, "kb.pdf"), width = 30, height = 10)
  
}

plotRawCoverage <- function(coverage_df, experiment_name, original_sample_order){
  
  # set sample order in facet
  coverage_df %<>%
    mutate(stage = factor(stage, levels = original_sample_order))
  
  # coverage plot
  ggplot(coverage_df, aes(x = pos, y = count)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    scale_fill_manual(values = c("grey", "black"), guide = FALSE) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + coverage_df$element_width[1])) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_grid(stage ~ .) +
    ggsave(filename = str_c(experiment_name, "_", coverage_df$element_name[1], "_minDistanceToGene_", coverage_df$min_distace[1] / 1000, "kb_sumCoverage.pdf"), width = 30, height = 10)
  
}
  
######################################################## READ DATA
# figure 6B and S5C
figure_6B_S5C_coverage_list <- readRDS(file = "figure_6B_S5C_coverage_list.RDS")

# figure S5D
figure_S5D_coverage_list <- readRDS(file = "figure_S5D_coverage_list.RDS")

# figure S5B
figure_S5B_coverage_list <- readRDS(file = "figure_S5B_coverage_list.RDS")

######################################################## MAIN CODE
### figure 6B and S5C 
# plot binned coverage
invisible(lapply(X = figure_6B_S5C_coverage_list, FUN = plotCummulativeFPKMBin, 
                 tile_width = 10000, 
                 experiment_name = "figure_6B_S5C", 
                 original_sample_order = c("GV", "2C", "2Ca"), 
                 plot_sample_order = c("2Ca", "2C", "GV")))
  

### figure S5D
# plot summed coverage
invisible(lapply(X = figure_S5D_coverage_list, FUN = plotRawCoverage, 
                 experiment_name = "figure_S5D", 
                 original_sample_order = c("GV", "1C", "2C", "2Ca")))


### figure S5B
# get data.frame from list
figure_S5B_coverage <- figure_S5B_coverage_list[[1]]

# plot summed coverage 
plotRawCoverage(coverage_df = figure_S5B_coverage, 
                experiment_name = "figure_S5B", 
                original_sample_order = c("s_Oo", "s_2C"))

# raw counts to FPM
figure_S5B_coverage_fpm <- 
  figure_S5B_coverage %>% 
  mutate(count = count / library_size) %>% 
  dplyr::select(-library_size)

# plot binned coverage
plotCummulativeFPKMBin(coverage_df = figure_S5B_coverage_fpm, 
                       tile_width = 10000, 
                       experiment_name = "figure_S5B", 
                       original_sample_order = c("s_Oo", "s_2C"), 
                       plot_sample_order = c("s_2C", "s_Oo"))