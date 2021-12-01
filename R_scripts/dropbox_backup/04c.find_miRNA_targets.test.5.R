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

# miR family info path
mir_info_path <- file.path(accessory_data_path, "targetScan", "miR_Family_Info.txt")

# all 3'UTR sequences 
aligned_3UTRs_path <- file.path(accessory_data_path, "aligned_UTR_sequences.mm.rn.bt.hs.RDS")

# miRNA annotated targets from targetScan path
mirna_targets_path <- file.path(accessory_data_path, "targetScan", "targetScan.mouse.7.1.20181109.conserved_family.conserved_targets.mouse.txt.gz")

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

# read miRNA gene targets from targetScan site
mirna_targets_table <- readr::read_delim(mirna_targets_path, 
                                         delim = "\t", 
                                         skip = 1, 
                                         col_names = c("mir_family", "gene_id", "gene_symbol", "transcript_id", "species_id", 
                                                       "UTR_start", "UTR_end", "MSA_start", "MSA_end", "seed_match", "PCT"))

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


##### 
### find miRNA family targets
# set top 10 miRNAs
mirna_top10 <- c("miR-128-3p", "miR-129-5p", "miR-140-3p.1", "miR-185-5p", "miR-186-5p", 
                 "miR-3064-5p/3085-3p", "miR-329-3p/362-3p", "miR-340-5p", "miR-495-3p", "miR-665-3p")

# loop through top 10 miRNA families
purrr::map(mirna_top10, function(unique_mir_family){
  
  # message
  cat("Plotting", unique_mir_family, "\n")
  
  # get targets from targetScan
  mirna_targets_targetScan <-
    mirna_targets_table %>%
    dplyr::filter(mir_family == unique_mir_family,
                  seed_match %in% c("7mer-m8", "8mer", "7mer-1a")) %>%
    arrange(gene_id)
  
  # get 3 patterns of one miRNA - 8mer, 7mer-m8 and 7mer-1a
  mir_patterns <- 
    mir_info %>% 
    filter(mir_family == unique_mir_family, 
           str_detect(mirbase_id, "^mmu-")) %>%
    dplyr::select(mature_sequence) %>% 
    dplyr::slice(1) %>% 
    dplyr::mutate(seed_7mer_m8 = str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T"), 
                  seed_8mer = str_c("T", str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T")), 
                  seed_7mer_1a = str_c("T", str_sub(mature_sequence, 2, 7) %>% str_replace_all(., "U", "T"))) %$%
    DNAStringSet(c(seed_7mer_m8, seed_8mer, seed_7mer_1a)) %>% 
    reverseComplement(.) %>% 
    magrittr::set_names(., c("7mer-m8", "8mer", "7mer-1a")) %>% 
    as.character(.) 
  
  # add gaps to pattern
  mir_patterns_gapped <- purrr::map(mir_patterns, function(pattern){
    
    # split, add gaps
    str_split(pattern, "", simplify = T) %>% 
      str_c(., collapse = "-*")
    
  }) 
  
  # match pattern on all animals 
  matches_list <- purrr::map(names(aligned_3UTRs), function(animal_name){
    
    # get animal 3'UTRs as character
    animal_3UTRs_char <- 
      aligned_3UTRs[[animal_name]] %>% 
      as.character(.) %>% 
      set_names(., names(aligned_3UTRs[[animal_name]]))
    
    # loop through 3 patterns
    mirna_targets <- purrr::map(names(mir_patterns_gapped), function(mir_pattern_name){
      
      # match
      mir_match <- stringi::stri_locate_all(animal_3UTRs_char, regex = mir_patterns_gapped[[mir_pattern_name]])
      
      # as table
      mir_match_tb <- 
        lapply(1:length(mir_match), function(n) cbind(mir_match[[n]], n)) %>% 
        do.call(rbind, .) %>% 
        as.tibble(.) %>% 
        dplyr::filter_at(vars("start", "end"), any_vars(!is.na(.))) %>% 
        dplyr::mutate(gene_id = names(animal_3UTRs_char)[n] %>% str_remove(., "\\|.*"),
                      pattern = mir_pattern_name) %>% 
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
  
  # get patterns which exist in aligned 3'UTR of all 4 animals
  mirna_targets_tb <- purrr::reduce(matches_list, dplyr::inner_join, by = c("gene_id", "start", "end", "pattern"))
  
  # add manual and targetScan miRNA targets to list
  mirna_targets_list <-
    list(mirna_targets_tb, mirna_targets_targetScan) %>% 
    set_names(., c("manual", "targetScan"))
  
  # plot manual and targetScan targets separately
  purrr::map(names(mirna_targets_list), function(target){
    
    # count targets, join with binned data
    mirna_targets <- 
      mirna_targets_list[[target]] %>% 
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
    heatmap_plot <- 
      ggplot(mirna_targets_smooth, aes(bin_absolute, bin_relative)) + 
      geom_tile(aes(fill = n)) + 
      scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"),
                           values = scales::rescale(c(-1, 0, 1))) +
      theme_tufte() 
    
    # save plot
    ggsave(plot = heatmap_plot, filename = file.path(outpath, str_c("grid.heatmap.Su.mas5", tissue, 
                                                                    "61.bins.targetScan", str_replace_all(unique_mir_family, "/", "-"), 
                                                                    "smooth", k, "conserved", target, "png", sep = ".")))
    
    # return
    return(str_c(unique_mir_family, target, sep = "_"))
    
  })
  
  # return  
  return(unique_mir_family)
  
})




