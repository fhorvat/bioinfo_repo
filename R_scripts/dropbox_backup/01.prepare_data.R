### INFO: 
### DATE: Sun Nov 04 19:08:17 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "gtfToGRanges.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# all 3'UTR sequences
all_3UTRs_path <- file.path(inpath, "targetScan", "targetScan.mouse.7.1.20181112.UTR_sequences.txt.gz")

# taxonomy ID path
tax_id_path <- file.path(inpath, "ncbi_taxonomy_names.clean.txt")

# miR family info path
mir_info_path <- file.path(inpath, "targetScan", "miR_Family_Info.txt")

######################################################## READ DATA
# read all aligned 3'UTRs from targetScan
all_3UTRs <- fread(all_3UTRs_path, col.names = c("refseq_id", "gene_id", "gene_symbol", "tax_id", "UTR_seq"), colClasses = "character")
all_3UTRs[tax_id == "9478", tax_id := "1868482"]
setkey(all_3UTRs, tax_id)

# read tax ID
tax_id <- fread(tax_id_path, col.names = c("tax_id", "tax_name", "unique_name", "name_class"), quote = "", colClasses = "character")
setkey(tax_id, tax_id)

# read miR info
mir_info <- readr::read_delim(mir_info_path, 
                              delim = "\t", 
                              skip = 1, 
                              col_names = c("mir_family", "seed_m8", "species_ID", "mirbase_id", "mature_sequence",
                                            "family_conservation", "mir_base_accession"))

######################################################## MAIN CODE
# #### DISTANCE GRID ####
# # set dimension
# dimension <- 50
# 
# # prepare full grid, save
# as.tibble(expand.grid(1:dimension, 1:dimension)) %>% 
#   set_colnames(c("bin_absolute", "bin_relative")) %>% 
#   dplyr::mutate(bin_id = str_c(bin_absolute, bin_relative, sep = ".")) %>% 
#   column_to_rownames("bin_id") %>%
#   dist(., method = "euclidean") %>% 
#   as.matrix(.) %>%
#   as.data.table(., keep.rownames = "bin_id") %>% 
#   melt(., id.vars = "bin_id", variable.name = "bin_id_dist", value.name = "dist") %T>%
#   saveRDS(., file = file.path(outpath, str_c("distance_grid.Euclidean.", dimension, "x", dimension, ".dt.RDS")))


#### ALIGNED 3' UTRs ####
## filter aligned 3'UTRs from targetScan (get only mouse, rat, cow and human) and save as fasta and in table
# clean tax_id table
tax_id %<>%
  .[tax_id %in% unique(all_3UTRs$tax_id) & name_class == "scientific name", .(tax_id, tax_name)] %>%
  .[, tax_name := str_replace_all(tax_name, " ", "_")] %>%
  .[]

# add taxon name to sequence table
all_3UTRs[tax_id, tax_name := i.tax_name, on = "tax_id"]

# get only human, rat, mouse and cow UTRs
filtered_3UTRs <- all_3UTRs[tax_name %in% c("Bos_taurus", "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus")]

# add tax name to gene_id, replace U with T
filtered_3UTRs[, `:=`(gene_id = str_c(gene_id, "|", tax_name),
                      UTR_seq = str_replace_all(UTR_seq, "U|u", "T"))]

# convert to DNAStringSet, save
filtered_3UTRs <-
  DNAStringSet(filtered_3UTRs$UTR_seq) %>%
  magrittr::set_names(., filtered_3UTRs$gene_id)

# # get list  DNAStringSets
# aligned_3UTRs <-
#   list("Mus_musculus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Mus_musculus")],
#        "Rattus_norvegicus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Rattus_norvegicus")],
#        "Bos_taurus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Bos_taurus")],
#        "Homo_sapiens" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Homo_sapiens")])

aligned_3UTRs <-
  list("Mus_musculus" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Mus_musculus")],
       "Homo_sapiens" = filtered_3UTRs[str_detect(names(filtered_3UTRs), "Homo_sapiens")])


#### SEED PATTERNS ####
# set vector with unique miR family IDs
# mirbase_ids <- c("mmu-let-7a-5p", "mmu-miR-290a-5p", "mmu-miR-122-5p", "hsa-miR-1-3p", "mmu-miR-9-5p")
mirbase_ids <- c("mmu-let-7a-5p", "mmu-miR-290a-5p", "mmu-miR-191-5p", "mmu-miR-134-5p", "mmu-miR-130a-3p", "mmu-miR-200a-3p")

# find targets of all miRNAs
mirna_targets <- purrr::map(mirbase_ids, function(unique_mirbase_id){
  
  # get 6mers and 7mers patterns of miRNA
  seed_patterns <- 
    mir_info %>% 
    filter(str_detect(mirbase_id, unique_mirbase_id)) %>%
    dplyr::select(mature_sequence) %>%
    dplyr::slice(1) %>% 
    dplyr::mutate(seed_6mer = str_sub(mature_sequence, 2, 7) %>% str_replace_all(., "U", "T"), 
                  seed_7mer = str_sub(mature_sequence, 2, 8) %>% str_replace_all(., "U", "T")) %$%
    DNAStringSet(c(seed_6mer, seed_7mer)) %>% 
    reverseComplement(.) %>% 
    as.character(.) %>% 
    magrittr::set_names(c("6mer", "7mer"))
  
  # add gaps to pattern
  seed_patterns_gapped <- purrr::map(seed_patterns, function(pattern){
    
    # split, add gaps
    str_split(pattern, "", simplify = T) %>% 
      str_c(., collapse = "-*")
    
  }) 
  
  # match pattern on all animals 
  seed_targets_list <- purrr::map(names(aligned_3UTRs), function(animal_name){
    
    # get animal 3'UTRs as character
    animal_3UTRs_char <- 
      aligned_3UTRs[[animal_name]] %>% 
      as.character(.) %>% 
      set_names(., names(aligned_3UTRs[[animal_name]]))
    
    # loop through seed patterns
    seed_targets <- purrr::map(names(seed_patterns_gapped), function(seed_pattern_name){
      
      # match
      seed_match <- stringi::stri_locate_all(animal_3UTRs_char, regex = seed_patterns_gapped[[seed_pattern_name]])
      
      # as table
      seed_match_tb <- 
        lapply(1:length(seed_match), function(n) cbind(seed_match[[n]], n)) %>% 
        do.call(rbind, .) %>% 
        as.tibble(.) %>% 
        dplyr::filter_at(vars("start", "end"), any_vars(!is.na(.))) %>% 
        dplyr::mutate(gene_id = names(animal_3UTRs_char)[n] %>% str_remove(., "\\|.*"),
                      pattern = seed_pattern_name) %>% 
        dplyr::select(-n)
      
    }) %>% 
      bind_rows(.) %>% 
      arrange(gene_id, start)
    
  })
  
  # get only targets which exist in aligned 3'UTR of all 4 animals ( = conserved targets), count targets
  seed_targets_tb <- 
    purrr::reduce(seed_targets_list, dplyr::inner_join, by = c("gene_id", "start", "end", "pattern")) %>% 
    dplyr::mutate(gene_id = str_remove(gene_id, "\\..*")) %>% 
    dplyr::count(gene_id, pattern, sort = T)
  
}) %>% 
  magrittr::set_names(., mirbase_ids) %T>%
  saveRDS(., file.path(outpath, "miRNA.mm_hsa.conserved_targets.mouse.6mers_7mers.RDS"))
