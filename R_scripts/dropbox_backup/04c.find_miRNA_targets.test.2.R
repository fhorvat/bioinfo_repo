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

# mouse 3' UTR sequences path
mouse_3UTRs_path <- file.path(accessory_data_path, "mm10.3pUTRs.DNAStringSet.RDS")

# # multiple alignment of 3' UTRs between mouse, rat, dog and human 
# aligned_3pUTRs_path <- file.path(accessory_data_path, "aligned_3pUTRs.mm10.rn5_hg19_canFam3.DNAStringSet.list.RDS")

# miR family info path
mir_info_path <- file.path(accessory_data_path, "targetScan", "miR_Family_Info.txt")

# mouse 3' UTR sequences targetScan
mouse_3UTRs_path_targetScan <- file.path(accessory_data_path, "targetScan", "targetScan.mouse.7.1.20181112.UTR_sequences.mouse.txt")

######################################################## READ DATA
# read distance data.table
dist_dt <- readRDS(distance_grid_path)

# read miRNA gene targets from targetScan site
mirna_targetScan <- readr::read_delim(mirna_targetScan_path, 
                                      delim = "\t", 
                                      skip = 1, 
                                      col_names = c("mir_family", "gene_id", "gene_symbol", "transcript_id", "species_id", 
                                                    "UTR_start", "UTR_end", "MSA_start", "MSA_end", "seed_match", "PCT"))

# read mouse 3' UTR sequences
mouse_3UTRs_seq <- readRDS(mouse_3UTRs_path)

# # read multiple alignment of 3' UTRs between mouse, rat, dog and human  
# aligned_3pUTRs <- readRDS(aligned_3pUTRs_path)

# read miR info
mir_info <- readr::read_delim(mir_info_path, 
                              delim = "\t", 
                              skip = 1, 
                              col_names = c("mir_family", "seed_m8", "species_ID", "mirbase_id", "mature_sequence",
                                            "family_conservation", "mir_base_accession"))

# read aligned UTRs from targetScan
mouse_3UTRs_seq_targetScan <- readr::read_delim(mouse_3UTRs_path_targetScan, 
                                                delim = "\t", 
                                                skip = 1, 
                                                col_names = c("refseq_id", "gene_id", "gene_symbol", "species_ID", "UTR_seq"))

######################################################## MAIN CODE
#####
## load grid
# set tissue (liver, skeletalmuscle) and miRNA to plot (skeletalmuscle = miR-133abc and miR-1ab/206/613, liver = miR-122)
tissue <- "skeletalmuscle"
mirna_family_id <- "miR-133-3p"

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
# ## get names of aligned UTRs
# aligned_3pUTRs_gene_id <- purrr::map(aligned_3pUTRs, function(aligned_3pUTR_animal){
#   
#   aligned_gene_id <- 
#     aligned_3pUTR_animal %>% 
#     names(.) %>% 
#     str_remove_all(., "^chr.*\\||\\.\\d+$")
#   
# })

##### 
## find all miRNA family targets
# set unique miR family
unique_mir_family <- "miR-133-3p"

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

### manual UTRs
# look for all 3 patterns in UTRs
mirna_targets_manual <- purrr::map(names(mir_pattern), function(mir_seed_name){
  
  # find one pattern
  Biostrings::vmatchPattern(pattern = mir_pattern[[mir_seed_name]], 
                            subject = mouse_3UTRs_seq, 
                            max.mismatch = 0) %>% 
    unlist(.) %>% 
    as.tibble(.) %>% 
    dplyr::select(-width) %>% 
    dplyr::mutate(seed_match = mir_seed_name)
  
}) %>% 
  bind_rows(.) %>% 
  arrange(names, start) %>% 
  dplyr::group_by(names, start) %>% 
  top_n(1, seed_match) %>% 
  dplyr::group_by(names, end) %>% 
  top_n(1, seed_match) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(gene_id = str_remove(names, "\\..*"))


### targets from targetScan
# get targets from targetScan
mirna_targets_targetScan <- 
  mirna_targetScan %>% 
  dplyr::filter(mir_family == unique_mir_family, 
                seed_match %in% c("7mer-m8", "8mer", "7mer-1a")) %>% 
  dplyr::filter((seed_match == "8mer" & PCT >= 0.6) | (seed_match == "7mer-m8" & PCT >= 1.8) | (seed_match == "7mer-1a" & PCT >= 2.5)) %>%
  arrange(gene_id)


### plot
# count targets
mirna_targets <- 
  mirna_targets_targetScan %>% 
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
ggplot(mirna_targets_smooth, aes(bin_absolute, bin_relative)) + 
  geom_tile(aes(fill = n)) + 
  scale_fill_gradientn(colours = c("darkblue", "blue", "green", "yellow", "orange", "red", "darkred"),
                       values = scales::rescale(c(-1, 0, 1))) +
  theme_tufte() +
  ggsave(filename = file.path(outpath,  str_c("grid.heatmap.Su.mas5", tissue, 
                                              "61.bins.targetScan", unique_mir_family, 
                                              "smooth", k, "conserved", "targetScan", "png", sep = ".")))


