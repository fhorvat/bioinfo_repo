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

# distance grid
distance_grid_path <- file.path(accessory_data_path, "distance_grid.Euclidean.50x50.dt.RDS")

# all 3'UTR sequences 
aligned_3UTRs_path <- file.path(accessory_data_path, "aligned_UTR_sequences.mm.rn.bt.hs.RDS")

# miR family info path
mir_info_path <- file.path(accessory_data_path, "targetScan", "miR_Family_Info.txt")

######################################################## READ DATA
# read distance data.table
dist_dt <- readRDS(distance_grid_path)

# read aligned 3'UTRs sequences in mouse, rat, cow and human
aligned_3UTRs <- readRDS(aligned_3UTRs_path)

# read miR info
mir_info <- readr::read_delim(mir_info_path, 
                              delim = "\t", 
                              skip = 1, 
                              col_names = c("mir_family", "seed_m8", "species_ID", "mirbase_id", "mature_sequence",
                                            "family_conservation", "mir_base_accession"))

######################################################## MAIN CODE
#####
## load grid and clean data
# set tissue (oocyte, liver, skeletalmuscle) and miRNA to plot
tissue <- "oocyte"

# grid path
grid_path <- file.path(inpath, str_c("grid.Su.mas5.", tissue, ".50.bins.csv"))

# read and clean grid
grid_df <- 
  readr::read_csv(grid_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))


##### 
### find septamer targets
# set unique miR family
unique_mirbase_id <- "mmu-let-7a-5p"

# get 3 patterns of one miRNA - 8mer, 7mer-m8 and 7mer-1a
septamer_patterns <- 
  mir_info %>% 
  filter(str_detect(mirbase_id, unique_mirbase_id)) %>%
  dplyr::select(mature_sequence) %>%
  dplyr::slice(1) %>% 
  dplyr::mutate(seed_7mer_m8 = str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T"), 
                seed_8mer = str_c("T", str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T")), 
                seed_7mer_1a = str_c("T", str_sub(mature_sequence, 2, 7) %>% str_replace_all(., "U", "T"))) %$%
  DNAStringSet(c(seed_7mer_m8, seed_8mer, seed_7mer_1a)) %>% 
  reverseComplement(.) %>% 
  as.character(.) %>% 
  magrittr::set_names(c("7mer-m8", "8mer", "7mer-1a"))

# get 3 patterns of one miRNA - 8mer, 7mer-m8 and 7mer-1a
# septamer_pattern <- "CTACCTC"
# septamer_patterns <-
#   tibble(seed_7mer_m8 = septamer_pattern) %>%
#   mutate(seed_7mer_1a = str_c(str_sub(seed_7mer_m8, 2, 7), "A")) %$%
#   c(seed_7mer_m8, seed_7mer_1a) %>%
#   magrittr::set_names(., c("7mer-m8", "7mer-1a"))

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
  dplyr::count(gene_id, sort = T) %>%
  dplyr::inner_join(., grid_df, by = "gene_id") %>%
  dplyr::group_by(bin_absolute, bin_relative) %>%
  dplyr::summarise(n = sum(n)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = ".")) %>% 
  as.data.table(., key = bin_id) %>% 
  .[order(bin_absolute, bin_relative), .(n, bin_id)] %>% 
  .[dist_dt, on = c("bin_id" = "bin_id_dist")] %>% 
  .[is.na(n), n := 0] %>% 
  .[]

# set smoothing factor k
k <- 20

# smooth 
septamer_targets_smooth <- 
  copy(septamer_targets_tb) %>% 
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
                                                                unique_mirbase_id, "smooth", k, "conserved", "png", 
                                                                sep = ".")))

# plot hits
heatmap_plot <- 
  ggplot(septamer_targets_smooth, aes(bin_absolute, bin_relative)) + 
  geom_point(aes(size = n)) + 
  # scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"),
  #                      values = scales::rescale(c(-1, 0, 1))) +
  theme_tufte() 

# save plot
ggsave(plot = heatmap_plot, filename = file.path(outpath, str_c("grid.heatmap.Su.mas5", tissue, "50.bins", 
                                                                unique_mirbase_id, "smooth", k, "conserved.dot", "png", 
                                                                sep = ".")))


