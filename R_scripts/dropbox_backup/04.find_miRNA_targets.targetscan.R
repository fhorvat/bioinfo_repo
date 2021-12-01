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

# miRNA gene targets from targetScan site path
mirna_targetScan_path <- file.path(accessory_data_path, "targetScan", "targetScan.mouse.7.1.20181109.conserved_family.conserved_and_unconserved_targets.mouse.txt.gz")

######################################################## READ DATA
# read distance data.table
dist_dt <- readRDS(distance_grid_path)

# read miRNA gene targets from targetScan site
mirna_targetScan <- readr::read_delim(mirna_targetScan_path, 
                                      delim = "\t", 
                                      col_names = c("mir_family", "gene_id", "gene_symbol", "transcript_id", "species_id", 
                                                    "UTR_start", "UTR_end", "MSA_start", "MSA_end", "seed_match", "PCT"), 
                                      skip = 1)

######################################################## MAIN CODE
#####
## load grid
# set tissue (liver, skeletalmuscle) and miRNA to plot (skeletalmuscle = miR-133abc and miR-1ab/206/613, liver = miR-122)
tissue <- "skeletalmuscle"

# grid path
grid_path <- file.path(inpath, str_c("grid.Su.mas5.", tissue, ".61.bins.csv"))

# read and clean grid
grid_df <- 
  readr::read_csv(grid_path) %>% 
  dplyr::select(gene_id, starts_with("bin"))


#####
## clean targetScan
mirna_targetScan %<>%
  dplyr::mutate(gene_id = str_remove(gene_id, "\\..*$"))


##### 
## create and plot grid
# get top 10 miRNA families by number of targets
mirna_top10 <- 
  mirna_targetScan %>% 
  count(mir_family) %>% 
  top_n(10, n) %$% 
  mir_family

# set miRNA family ID
# mirna_family_id <- "miR-133-3p"

## loop through top 10 miRNA families
purrr::map(mirna_top10, function(mirna_family_id){
  
  # find 3' UTR with miRNA targets, get only conserved (PCT: 8mer >= 0.6; 7mer-m8 >= 1.8; 7mer-1A >= 2.5)
  mirna_targets <- 
    mirna_targetScan %>% 
    dplyr::filter(mir_family == mirna_family_id, 
                  seed_match %in% c("7mer-m8", "8mer", "7mer-1a")) %>% 
    dplyr::filter((seed_match == "8mer" & PCT >= 0.6) | (seed_match == "7mer-m8" & PCT >= 1.8) | (seed_match == "7mer-1a" & PCT >= 2.5)) %>%
    dplyr::count(gene_id, sort = T) %>% 
    dplyr::inner_join(., grid_df, by = "gene_id") %>% 
    dplyr::group_by(bin_absolute, bin_relative) %>% 
    dplyr::summarise(n = sum(n)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = "."))
  
  # plot heatmap
  mirna_targets_dt <-
    mirna_targets %>% 
    as.data.table(., key = bin_id) %>% 
    .[order(bin_absolute, bin_relative), .(n, bin_id)] %>% 
    .[dist_dt, on = c("bin_id" = "bin_id_dist")] %>% 
    .[is.na(n), n := 0] %>% 
    .[]
  
  # set smoothing factor k
  k <- 25
  
  # smooth 
  mirna_targets_smooth <- 
    copy(mirna_targets_dt) %>% 
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
                                                "61.bins.targetScan", str_replace_all(mirna_family_id, "/", "-"), 
                                                "smooth", k, "conserved", "png", sep = ".")))
  
})




