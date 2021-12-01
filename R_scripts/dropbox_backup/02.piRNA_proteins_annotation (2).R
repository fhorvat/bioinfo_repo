### INFO: 
### DATE: Sun Oct 06 17:21:37 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_proteins_expression/results")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genomes_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_proteins_expression/genomes"

# gene info path
genes_info_paths <- list.files(path = genomes_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames..*geneInfo.csv$"), full.names = T, recursive = T)
genes_info_paths <- genes_info_paths[!str_detect(genes_info_paths, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.geneInfo.csv")]

######################################################## READ DATA
# read genes info
genes_info <- 
  purrr::map(genes_info_paths, readr::read_csv) %>% 
  set_names(., str_extract(genes_info_paths, "cow|mouse|rat|golden_hamster"))

######################################################## MAIN CODE
# Aravin et al. (10.1371/journal.pgen.1000764): 
# MILI-TDRD1 - The first type of granules, pi-bodies, contains the MILI-TDRD1 module of the piRNA pathway
# MIWI2-TDRD9-MAEL - The second type of granules, piP-bodies, harbors the MIWI2-TDRD9-MAEL module of the piRNA 
# MVH - mouse VASA homolog (MVH) protein, an RNA helicase
# 
# Ernst et al. (10.1038/s41467-017-01049-7): 
# Maelstrom (MAEL) - piRNA precursor molecules are then exported into the cytoplasm, which might be facilitated by Maelstrom (MAEL)
# PLD6 (MITOPLD) - the endonuclease that generates the 5' end of primary piRNAs = Zuccchini
# Tdrkh - The process is facilitated by Tdrkh, a Tudor and KH domain-containing protein, whose absence causes the accumulation of 31–37 nucleotide (nt) long intermediates
# HENMT1 - Mature 3' ends are then 2'-O methylated by the RNA methyltransferase HENMT1
# RNF17 - The disengagement of MIWI from the ping-pong cycle was shown to be due to active repression of the ping-pong cycle by RNF17, a Tudor family member, in meiotic cells
# A-MYB - Transcription of these clusters is facilitated by the ancestral transcription factor A-MYB
# DDX4 - Vasa (also known as Ddx4) is an ATP-dependent RNA helicase

# # Pandey et al. (10.1073/pnas.1316316110)
# TDRD12 - Tudor domain containing 12 (TDRD12) is essential for secondary PIWI interacting RNA biogenesis in mice
# 
# Vagin et al. (10.1101/gad.1814809): 
# Piwi proteins were found in complex with PRMT5/WDR77, an enzyme that dimethylates arginine residues
# 
# Chen et al. (10.1073/pnas.0911640106): 
# Tdrd1, Tdrd4/RNF17, and Tdrd6, Tdrkh/Tdrd2
# In this analysis, we observed several previously known Piwi-interacting proteins or piRNA pathway components such as Ddx4, Kif17, Fxr1, and Mov10l1

# Siomi et al. (10.1101/gad.1899210)
# PIWI proteins have been shown recently to contain symmetrical dimethyl arginines (sDMAs), 
# and this modification is mediated by the methyltransferase PRMT5 (also known as Dart5 or Capsuleen)

### piRNA-pathway associated proteins
# create table
pirna_genes_tb <- 
  tibble(gene_name = c("Piwil2", "Piwil1", "Piwil4", "Piwil3",
                       "Mov10l1", "Ddx4", "Mael", 
                       "Pld6", "Henmt1",
                       "Tdrd1", "Tdrkh", "Rnf17", "Tdrd6", "Tdrd9", "Tdrd12", 
                       "Prmt5", "Wdr77", "Kif17"),
         other_name = c("MILI", "MIWI", "MIWI2", "MIWI3",
                        "MOV10L1", "DDX4/VASE", "MAEL/Maelstorm", 
                        "PLD6/MITOPOLD/Zucchini", "HENMT1", 
                        "TDRD1", "TDRD2/TDRKH", "TDRD4/RNF17", "TDRD6", "TDRD9", "TDRD12", 
                        "PRMT5/Dart5/Capsuleen", "WDR77", "KIF17"))

# filter gene names
gene_names_list <- 
  purrr::map(names(genes_info), function(animal){
    
    # filter
    gene_names_tb <- genes_info[[animal]]

    # if cow then uppercase gene names
    if(animal == "cow"){
      
      gene_names_tb %<>%
        dplyr::filter(gene_name %in% toupper(pirna_genes_tb$gene_name)) %>% 
        dplyr::mutate(gene_name = stringr::str_to_title(gene_name))
      
    }else{
      
      gene_names_tb %<>%
        dplyr::filter(gene_name %in% pirna_genes_tb$gene_name)
      
    }
    
    # tidy
    gene_names_tb %<>%
      dplyr::select(gene_name, gene_id) %>% 
      set_names(., c("gene_name", str_c("gene_id", animal, sep = ".")))
    
    # return
    return(gene_names_tb)
    
  }) %>% 
  set_names(names(genes_info))

# join
gene_names_tb <- 
  purrr::reduce(gene_names_list, full_join, by = "gene_name") %>% 
  dplyr::select(gene_name, gene_id.mouse, gene_id.rat, gene_id.golden_hamster, gene_id.cow)

# save
readr::write_csv(gene_names_tb, path = file.path(outpath, "piRNA_proteins.gene_ids.csv"))


### coordinates and save
# add coordinates
gene_coordinates <- 
  readr::read_csv(file = file.path(outpath, "piRNA_proteins.gene_ids.with_notes.csv")) %>% 
  dplyr::select(-note) %>% 
  pivot_longer(., cols = -gene_name, names_to = "animal", values_to = "gene_id", names_prefix = "gene_id.") %>% 
  dplyr::left_join(., genes_info %>% bind_rows(.) %>% dplyr::select(gene_id, seqnames, start, end, strand), by = "gene_id") %>% 
  tidyr::unite(col = "coordinates", c("seqnames", "start", "end"), sep = " ")

# save
readr::write_csv(gene_coordinates, path = file.path(outpath, "piRNA_proteins.coordinates.csv"))

