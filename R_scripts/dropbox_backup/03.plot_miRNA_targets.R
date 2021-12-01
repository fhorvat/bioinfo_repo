### INFO: 
### DATE: Sun Nov 04 19:08:17 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/biogps")

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

library(ggthemes)
library(Biostrings)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# plots denisty heatmap on binned grid
plotHeatmap <- function(df, scale = F, min = 0, max = 10, save = F){
  
  # get tissue and miRNA
  tissue <- unique(df$tissue)
  mirna <- unique(df$mirna)
  
  # plot
  g <- ggplot(df, aes(bin_absolute, bin_relative)) + 
    geom_tile(aes(fill = n)) + 
    theme_tufte() + 
    labs(title = str_c(tissue, " ", mirna)) + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())
  
  # scale
  if(scale){
    
    # add scale
    g <- 
      g +
      scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"), 
                           limits = c(min, max))
    
    # save
    if(save){
      g + suppressMessages(ggsave(filename = file.path(outpath, "individual_plots", "scaled", 
                                                       str_c("grid.biogps.gcrma", "scaled", 
                                                             tissue, "50.bins", mirna, kmer, "smooth", k, 
                                                             "conserved", "png", sep = "."))))
    }
    
  }else{
    
    # add scale without scaling
    g <- 
      g +
      scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"))
    
    # save
    if(save){
      g + suppressMessages(ggsave(filename = file.path(outpath, "individual_plots", "not_scaled", 
                                                       str_c("grid.biogps.gcrma", "not_scaled", 
                                                             tissue, "50.bins", mirna, kmer, "smooth", k, 
                                                             "conserved", "png", sep = "."))))
    }
    
  }
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# accessory data path
accessory_data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data"

# distance grid
distance_grid_path <- file.path(accessory_data_path, "distance_grid.Euclidean.50x50.dt.RDS")

# binned grid path
binned_grid_path <- list.files(inpath, "grid\\.biogps\\.gcrma.*RDS")

# miRNA conserved 6mers and 7mers targets in 3'UTRs of mouse 
mirna_targets_path <- file.path(accessory_data_path, "miRNA.mm_rn_bt_hsa.conserved_targets.mouse.6mers_7mers.RDS")

######################################################## READ DATA
# read distance data.table
dist_dt <- readRDS(distance_grid_path)

# read and clean binned grids
binned_grid_list <- 
  purrr::map(binned_grid_path, function(grid){
    
    # read binned grid
    readRDS(grid) %>% 
      dplyr::select(gene_id, starts_with("bin"))
    
  }) %>% 
  magrittr::set_names(., str_remove_all(binned_grid_path, "grid\\.biogps\\.gcrma\\.|\\.[0-9]{2}.*$"))

# read miRNA targets
mirna_targets <- readRDS(mirna_targets_path) 

######################################################## MAIN CODE
# set tissue and miRNAs to check targeting of specific miRNA
tissue_list_unique <- c("oocyte", "fertilizedegg", "blastocysts", "embryonic_stem_no_feeder", "embryonic_stem_feeder_layer", "liver", "skeletalmuscle")
mirna_list_unique <- c("mmu-let-7a-5p", "mmu-miR-290a-5p", "mmu-miR-191-5p", "mmu-miR-134-5p", "mmu-miR-130a-3p", "mmu-miR-200a-3p")

# replicate tissue list
tissue_list <- rep(tissue_list_unique, length(mirna_list_unique))

# set miRNA-tissue relations
tissue_mirna_list <- 
  rep(mirna_list_unique, each = length(tissue_list_unique)) %>% 
  magrittr::set_names(., tissue_list)

### plot separately for 6mers and 7mers
# set smoothing factor k
k <- 5

# get density of 6mers/7mers target in different tissue
kmer_binned_density <- purrr::map(c("6mer", "7mer"), function(kmer){
  
  # plot kmer in all tissues
  seed_targets_smooth <- purrr::map(1:length(tissue_list), function(n){
    
    # get tissue
    tissue <- names(tissue_mirna_list[n])
    
    # get unique target miRNA
    tissue_mirna <- tissue_mirna_list[n]
    
    # message
    cat(tissue, tissue_mirna, "\n")
    
    # get binned targets, smooth
    seed_targets_smooth <- 
      mirna_targets[[tissue_mirna]] %>%
      dplyr::filter(pattern == kmer) %>% 
      dplyr::inner_join(., binned_grid_list[[tissue]], by = "gene_id") %>%
      dplyr::group_by(bin_absolute, bin_relative) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = ".")) %>% 
      as.data.table(., key = bin_id) %>% 
      .[order(bin_absolute, bin_relative), .(n, bin_id)] %>% 
      .[dist_dt, on = c("bin_id" = "bin_id_dist")] %>% 
      .[is.na(n), n := 0] %>% 
      .[, dist := (1 / (dist ^ 2 + k))] %>% 
      .[, n := (n * dist)] %>% 
      .[, .(n = sum(n)), by = i.bin_id] %>% 
      .[, c("bin_absolute", "bin_relative") := lapply(tstrsplit(i.bin_id, ".", fixed = TRUE), as.integer)] %>% 
      .[order(bin_absolute, bin_relative)] %>% 
      .[, `:=`(tissue = tissue, mirna = tissue_mirna)] %>% 
      .[]
    
  })
  
}) %>% 
  magrittr::set_names(., c("6mer", "7mer"))


### plot 
# plot tissue/miRNA pairs for 6mers/7mers
purrr::map(c("6mer", "7mer"), function(kmer){
  
  # get all densities in one table
  kmer_density_all <- 
    bind_rows(kmer_binned_density[[kmer]]) %>% 
    as.tibble(.) %>% 
    mutate(tissue_mirna = str_c(tissue, ".", mirna),
           tissue_mirna = factor(tissue_mirna, levels = str_c(names(tissue_mirna_list), tissue_mirna_list, sep = ".")))
  
  # get range of colors for plots
  min_range <- min(kmer_density_all$n)
  max_range <- max(kmer_density_all$n)
  
  # get plots in tibble
  kmer_density_all %<>% 
    group_by(tissue_mirna) %>% 
    nest(.) %>% 
    mutate(plots_not_scaled = map(data, plotHeatmap, scale = F, save = T), 
           plots_scaled = map(data, plotHeatmap, scale = T, save = T, min = min_range, max = max_range))
  
  # plot as panel
  gridExtra::arrangeGrob(grobs = kmer_density_all$plots_not_scaled, ncol = 7) %>% 
    ggsave(plot = ., filename = file.path(outpath, 
                                          str_c("grid.biogps.gcrma", "not_scaled", "50.bins", 
                                                kmer, "smooth", k, 
                                                "conserved", "png", sep = ".")), width = 40, height = 30)
  
  # plot as panel
  gridExtra::arrangeGrob(grobs = kmer_density_all$plots_scaled, ncol = 7) %>% 
    ggsave(plot = ., filename = file.path(outpath, 
                                          str_c("grid.biogps.gcrma", "scaled", "50.bins", 
                                                kmer, "smooth", k, 
                                                "conserved", "png", sep = ".")), width = 40, height = 30)
  
  # return
  return(kmer)
  
})
