### INFO: 
### DATE: Sun Nov 04 19:08:17 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133")

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

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# accessory data path
accessory_data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data"

# all 3'UTR sequences 
aligned_3UTRs_path <- file.path(accessory_data_path, "aligned_UTR_sequences.mm.rn.bt.hs.RDS")

######################################################## READ DATA
# read aligned 3'UTRs sequences in mouse, rat, cow and human
aligned_3UTRs <- readRDS(aligned_3UTRs_path)

######################################################## MAIN CODE
#####
## load grid and clean data
# set tissue (oocyte, liver, skeletalmuscle) and miRNA to plot
tissue <- "oocyte"


### binned grid
# grid path
grid_array_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133/grid.Su.mas5.oocyte.50.bins.filtered.csv")
grid_ngs_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/Su_2004_ProcNatlAcadSciUSA_GSE1133/grid.Fugaku.Encode.oocyte.50.bins.filtered.csv")

# read and clean array grid
grid_array <- 
  readr::read_csv(grid_array_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))

# read and clean NGS grid
grid_ngs <- 
  readr::read_csv(grid_ngs_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))

# add to list
grid_list <- list(array = grid_array, ngs = grid_ngs)


### distance grid
# distance grid paths
distance_array_path <- file.path(accessory_data_path, "distance_grid.Euclidean.50x50.dt.RDS")
distance_ngs_path <- file.path(accessory_data_path, "distance_grid.Euclidean.50x50.dt.RDS")

# read distance data.table
dist_array_dt <- readRDS(distance_array_path)
dist_ngs_dt <- readRDS(distance_ngs_path)

# add to list
distance_list <- list(array = dist_array_dt, ngs = dist_ngs_dt)
  
##### 
### find miRNA family targets
# set top 10 septamers
septamer_top10 <- c("miR-128-3p", "miR-129-5p", "miR-140-3p.1", "miR-185-5p", "miR-186-5p", 
                    "miR-3064-5p/3085-3p", "miR-329-3p/362-3p", "miR-340-5p", "miR-495-3p", "miR-665-3p")

# loop through top 10 miRNA families
purrr::map(septamer_top10, function(septamer_pattern){
  
  septamer_pattern <- "CACTGTG"
  
  # message
  cat("Plotting", septamer_pattern, "\n")
  
  # get 3 patterns of one miRNA - 8mer, 7mer-m8 and 7mer-1a
  septamer_patterns <- 
    tibble(seed_7mer_m8 = septamer_pattern) %>% 
    mutate(seed_7mer_1a = str_c(str_sub(seed_7mer_m8, 2, 7), "A")) %$%
    c(seed_7mer_m8, seed_7mer_1a) %>% 
    magrittr::set_names(., c("7mer-m8", "7mer-1a"))
  
  # add gaps to pattern
  septamer_patterns_gapped <- purrr::map(septamer_patterns, function(pattern){
    
    # split, add gaps
    str_split(pattern, "", simplify = T) %>% 
      str_c(., collapse = "-*")
    
  }) 
  
  # match pattern on all animals 
  septamer_targets_list <- purrr::map(names(aligned_3UTRs), function(animal_name){
    
    # get animal 3'UTRs as character
    animal_3UTRs_char <- 
      aligned_3UTRs[[animal_name]] %>% 
      as.character(.) %>% 
      set_names(., names(aligned_3UTRs[[animal_name]]))
    
    # loop through 3 patterns
    septamer_targets <- purrr::map(names(septamer_patterns_gapped), function(septamer_pattern_name){
      
      # match
      septamer_match <- stringi::stri_locate_all(animal_3UTRs_char, regex = septamer_patterns_gapped[[septamer_pattern_name]])
      
      # as table
      septamer_match_tb <- 
        lapply(1:length(septamer_match), function(n) cbind(septamer_match[[n]], n)) %>% 
        do.call(rbind, .) %>% 
        as.tibble(.) %>% 
        dplyr::filter_at(vars("start", "end"), any_vars(!is.na(.))) %>% 
        dplyr::mutate(gene_id = names(animal_3UTRs_char)[n] %>% str_remove(., "\\|.*"),
                      pattern = septamer_pattern_name) %>% 
        dplyr::select(-n)
      
    }) %>% 
      bind_rows(.) %>% 
      arrange(gene_id, start) %>% 
      dplyr::group_by(gene_id, start) %>% 
      top_n(1, pattern) %>% 
      dplyr::group_by(gene_id, end) %>% 
      top_n(1, pattern) %>% 
      dplyr::ungroup(.)
    
  })
  
  # get only targets which exist in aligned 3'UTR of all 4 animals (= conserved targets), count targets, join with binned data
  septamer_targets_tb <- 
    purrr::reduce(septamer_targets_list, dplyr::inner_join, by = c("gene_id", "start", "end", "pattern")) %>% 
    dplyr::mutate(gene_id = str_remove(gene_id, "\\..*")) %>% 
    dplyr::count(gene_id, sort = T)
  
  # plot on NGS and array grids
  purrr::map(names(grid_list), function(grid){
    
    # join with grid table
    septamer_targets_grid <- 
      septamer_targets_tb %>% 
      dplyr::inner_join(., grid_list[[grid]], by = "gene_id") %>%
      dplyr::group_by(bin_absolute, bin_relative) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = ".")) %>% 
      as.data.table(., key = bin_id) %>% 
      .[order(bin_absolute, bin_relative), .(n, bin_id)] %>% 
      .[distance_list[[grid]], on = c("bin_id" = "bin_id_dist")] %>% 
      .[is.na(n), n := 0] %>% 
      .[]
    
    # set smoothing factor k
    k <- 25
    
    # smooth 
    septamer_targets_smooth <- 
      copy(septamer_targets_grid) %>% 
      .[, dist := (1 / (dist ^ 2 + k))] %>% 
      .[, n := (n * dist)] %>% 
      .[, .(n = sum(n)), by = i.bin_id] %>% 
      .[, n := as.vector(scale(n))] %>% 
      .[, c("bin_absolute", "bin_relative") := lapply(tstrsplit(i.bin_id, ".", fixed = TRUE), as.integer)] %>% 
      .[order(bin_absolute, bin_relative)]
    
    # plot hits
    heatmap_plot <- 
      ggplot(septamer_targets_smooth, aes(bin_absolute, bin_relative)) + 
      geom_tile(aes(fill = n)) + 
      scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"),
                           values = scales::rescale(c(-1, 0, 1))) +
      theme_tufte() 
    
    # save plot
    ggsave(plot = heatmap_plot, filename = file.path(outpath, str_c("grid.heatmap.Su.mas5", tissue, "50.bins", 
                                                                    septamer_pattern, grid, "smooth", k, "conserved", "png", 
                                                                    sep = ".")))
    
  })
  
  
  # return
  return(septamer_pattern)
  
})





