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

library(ggthemes)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
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
distance_grid_path <- file.path(accessory_data_path, "distance_grid.Euclidean.61x61.dt.RDS")

# mouse 3' UTR sequences path
mouse_3UTRs_path <- file.path(accessory_data_path, "mm10.3pUTRs.DNAStringSet.RDS")

# mouse 3' UTR sequences targetScan
mouse_3UTRs_path_targetScan <- file.path(accessory_data_path, "targetScan", "targetScan.mouse.7.1.20181112.UTR_sequences.mouse.txt")

# miR family info path
mir_info_path <- file.path(accessory_data_path, "targetScan", "miR_Family_Info.txt")

######################################################## READ DATA
# read distance data.table
dist_dt <- readRDS(distance_grid_path)

# read mouse 3'UTR sequences from ENSEMBL
mouse_3UTRs_ensembl <- readRDS(mouse_3UTRs_path)

# read mouse 3'UTRs sequences from targetScan
mouse_3UTRs_targetScan <- readr::read_delim(mouse_3UTRs_path_targetScan, 
                                                   delim = "\t", 
                                                   skip = 1, 
                                                   col_names = c("refseq_id", "gene_id", "gene_symbol", "species_ID", "UTR_seq"))

# read miR info
mir_info <- readr::read_delim(mir_info_path, 
                              delim = "\t", 
                              skip = 1, 
                              col_names = c("mir_family", "seed_m8", "species_ID", "mirbase_id", "mature_sequence",
                                            "family_conservation", "mir_base_accession"))

######################################################## MAIN CODE
#####
## load grid and clean data
# set tissue (liver, skeletalmuscle) and miRNA to plot
tissue <- "skeletalmuscle"

# grid path
grid_path <- file.path(inpath, str_c("grid.Su.mas5.", tissue, ".61.bins.csv"))

# read and clean grid
grid_df <- 
  readr::read_csv(grid_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))

# clean mouse 3'UTRs sequences from targetScan
mouse_3UTRs_targetScan <-
  mouse_3UTRs_targetScan %$%
  UTR_seq %>%
  str_replace_all(., "U|u", "T") %>%
  DNAStringSet(.) %>%
  DECIPHER::RemoveGaps(.) %>%
  set_names(., mouse_3UTRs_targetScan$gene_id)


##### 
### find miRNA family targets
# set top 10 miRNAs
mirna_top10 <- c("miR-128-3p", "miR-129-5p", "miR-140-3p.1", "miR-185-5p", "miR-186-5p", 
                 "miR-3064-5p/3085-3p", "miR-329-3p/362-3p", "miR-340-5p", "miR-495-3p", "miR-665-3p")

# loop through top 10 miRNA families
purrr::map(mirna_top10, function(unique_mir_family){
  
  # get 3 patterns of one miRNA - 8mer, 7mer-m8 and 7mer-1a
  mir_pattern <- 
    mir_info %>% 
    filter(str_detect(mir_family, unique_mir_family), 
           str_detect(mirbase_id, "^mmu-")) %>%
    dplyr::select(mature_sequence) %>% 
    dplyr::slice(1) %>% 
    dplyr::mutate(seed_7mer_m8 = str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T"), 
                  seed_8mer = str_c("T", str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T")), 
                  seed_7mer_1a = str_c("T", str_sub(mature_sequence, 2, 7) %>% str_replace_all(., "U", "T"))) %$%
    DNAStringSet(c(seed_7mer_m8, seed_8mer, seed_7mer_1a)) %>% 
    reverseComplement(.) %>% 
    set_names(c("7mer-m8", "8mer", "7mer-1a"))
  
  # look for seed patterns in 3 UTRs downloaded from ENSEMBL
  mirna_targets_ensembl <- purrr::map(names(mir_pattern), function(mir_seed_name){
    
    # find one pattern
    Biostrings::vmatchPattern(pattern = mir_pattern[[mir_seed_name]], 
                              subject = mouse_3UTRs_ensembl, 
                              max.mismatch = 0) %>% 
      unlist(.) %>% 
      as.tibble(.) %>% 
      dplyr::select(-width) %>% 
      dplyr::mutate(seed_match = mir_seed_name)
    
  })
  
  # look for seed patterns in 3'UTRs downloaded from targetScan
  mirna_targets_targetScan <- purrr::map(names(mir_pattern), function(mir_seed_name){
    
    # find one pattern
    Biostrings::vmatchPattern(pattern = mir_pattern[[mir_seed_name]], 
                              subject = mouse_3UTRs_targetScan, 
                              max.mismatch = 0) %>% 
      unlist(.) %>% 
      as.tibble(.) %>% 
      dplyr::select(-width) %>% 
      dplyr::mutate(seed_match = mir_seed_name)
    
  })
  
  # put ensembl and targetScan to list
  mirna_targets_list <- 
    list(mirna_targets_ensembl, mirna_targets_targetScan) %>% 
    set_names(., c("ensembl", "targetScan"))
  
  ### plot miRNA targets for ENSEMBL and targetScan UTRs
  purrr::map(names(mirna_targets_list), function(mirna_targets_utr){
    
    # get all targets, join with binned genes
    mirna_targets <- 
      mirna_targets_list[[mirna_targets_utr]] %>% 
      bind_rows(.) %>% 
      arrange(names, start) %>% 
      dplyr::group_by(names, start) %>% 
      top_n(1, seed_match) %>% 
      dplyr::group_by(names, end) %>% 
      top_n(1, seed_match) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::mutate(gene_id = str_remove(names, "\\..*")) %>% 
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
    k <- 25
    
    # smooth 
    mirna_targets_smooth <- 
      copy(mirna_targets) %>% 
      .[, dist := (1 / (dist ^ 2 + k))] %>% 
      .[, n := (n * dist)] %>% 
      .[, .(n = sum(n)), by = i.bin_id] %>% 
      .[, n := as.vector(scale(n))] %>% 
      .[, c("bin_absolute", "bin_relative") := lapply(tstrsplit(i.bin_id, ".", fixed = TRUE), as.integer)] %>% 
      .[order(bin_absolute, bin_relative)]
    
    # plot hits
    ggplot(mirna_targets_smooth, aes(bin_absolute, bin_relative)) + 
      geom_tile(aes(fill = n)) + 
      scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"),
                           values = scales::rescale(c(-1, 0, 1))) +
      theme_tufte() +
      ggsave(filename = file.path(outpath,  str_c("grid.heatmap.Su.mas5", tissue, 
                                                  "61.bins.targetScan", str_replace_all(unique_mir_family, "/", "-"), 
                                                  "smooth", k, "all", "manual", mirna_targets_utr, "png", sep = ".")))
    
  })
  
})

