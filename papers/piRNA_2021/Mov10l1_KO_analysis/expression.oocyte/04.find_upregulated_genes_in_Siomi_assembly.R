### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/Mov10l1_KO_analysis/oocyte.RNA_seq")

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

library(GenomicRanges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# diff. exp. results path
results_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression"
results_path <- list.files(results_path, ".*significant_results\\.xlsx", recursive = T, full.names = T)


### Siomi's genome
# base path
siomi_path <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# RefSeq .gtf lifted to Siomi's 
refseq_gtf_path <- file.path(siomi_path, "annotation/Liftoff/MesAur1/RefSeq")
refseq_gtf_path <- list.files(refseq_gtf_path, ".*\\.liftoff\\.gff", full.names = T)

# Siomi's repeatMasker
rmsk_clean_path <- list.files(siomi_path, "rmsk\\..*\\.clean\\.fa\\.out\\.gz", full.names = T)
rmsk_joined_path <- list.files(siomi_path, "rmsk\\..*\\.joined_rmsk_id\\.fa\\.out\\.gz", full.names = T)


### full length LINE1/IAP insertions 
# full length path
fli_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements"
fli_paths <- list.files(fli_path, pattern = ".*\\.FLI_elements\\.bed", full.names = T, recursive = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read diff. exp. results
results_tb <- 
  openxlsx::read.xlsx(results_path, sheet = "Mov10l_KO_vs_Mov10l_WT") %>% 
  as_tibble(.)

# read RefSeq gtf
refseq_gtf_tb <- readr::read_delim(file = refseq_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))
refseq_gtf_gr <- rtracklayer::import.gff(refseq_gtf_path)

# read rmsk
rmsk_tb <- readr::read_delim(rmsk_joined_path, delim = "\t")

# read full length insertions
fli_list <- 
  purrr::map(fli_paths, rtracklayer::import.bed) %>% 
  set_names(., str_extract(fli_paths, "LINE1|IAP"))

######################################################## MAIN CODE
### clean files
# clean Siomi's gtf
gtf_genes <- 
  refseq_gtf_gr %>% 
  as_tibble(.) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::select(seqnames, start, end, strand, ID, Dbxref, Name, gene_biotype) %>% 
  dplyr::mutate(Dbxref = unlist(Dbxref), 
                Dbxref = str_remove(Dbxref, "GeneID:")) %>% 
  tidyr::unite(col = coordinates, seqnames, start, end, strand, sep = " ")

# get average length of full length insertions, separately for IAP and LINE1
fli_tb <- purrr::map(names(fli_list), function(name){
    
    # find mean insertion length
    fli_mean <- 
      fli_list[[name]] %>% 
      as_tibble(.) %>% 
      dplyr::summarize(mean_width = mean(width)) %>% 
      dplyr::mutate(element = name) %>% 
      dplyr::select(element, mean_width)
    
    # return
    return(fli_mean)
    
  }) %>% 
  dplyr::bind_rows(.)

# get GRanges repeatMasker, filter LINE1s and IAPs which are at max. 0.5 kb shorter than mean width of FLIs
fli_gr <- 
  rmsk_tb %>% 
  dplyr::filter((repFamily == "L1") | (str_detect(repName, "^IAP"))) %>% 
  dplyr::mutate(element = ifelse(repClass == "LINE", "LINE1", "IAP")) %>% 
  dplyr::mutate(width = end - start + 1) %>% 
  dplyr::left_join(., fli_tb, by = "element") %>% 
  dplyr::filter(width >= (mean_width - 500)) %>% 
  dplyr::mutate(strand = ifelse(strand %in% c("+", "-"), strand, "*")) %>% 
  dplyr::select(-c(width, mean_width, element)) %>% 
  dplyr::mutate(repClass = ifelse(repClass == "LINE", "LINE1", "IAP")) %>% 
  GRanges(.)


### join by gene name
# # clean genes
# gene_names <- 
#   gtf_genes %>% 
#   dplyr::select(gene_name = Name, coordinates.Siomi = coordinates)
# 
# # filter results
# results_up <- 
#   results_tb %>% 
#   dplyr::filter(log2FoldChange >= 0) %>% 
#   dplyr::left_join(., gene_names, by = "gene_name") %>% 
#   dplyr::select(gene_id, gene_name, coordinates.mesAur1 = coordinates, coordinates.Siomi,
#                 log2FoldChange, padj, Mov10l_WT.FPKM, Mov10l_KO.FPKM,
#                 gene_biotype, gene_description)

### join by NCBI accession 
# clean genes
gene_names <- 
  gtf_genes %>% 
  dplyr::select(Dbxref, coordinates.Siomi = coordinates)

# filter results
results_up <- 
  results_tb %>% 
  dplyr::filter(log2FoldChange >= 0) %>%
  dplyr::mutate(Dbxref = str_remove_all(gene_description, ".*\\[.*Acc:|\\]")) %>% 
  dplyr::left_join(., gene_names, by = "Dbxref") %>% 
  dplyr::select(gene_id, gene_name, Dbxref, coordinates.mesAur1 = coordinates, coordinates.Siomi,
                log2FoldChange, padj, Mov10l_WT.FPKM, Mov10l_KO.FPKM,
                gene_biotype, gene_description)


### find nearest full length insertion to each gene
# get coordinates in Siomi's 
results_gr <- 
  results_up %>% 
  dplyr::filter(!is.na(Dbxref)) %>% 
  tidyr::separate(coordinates.Siomi, into = c("seqnames", "start", "end", "strand"), sep = " ") %>% 
  GRanges(.)

# find nearest FLIs to each gene
nearest_fli <- GenomicRanges::distanceToNearest(results_gr, fli_gr, ignore.strand = T)

# create results table
results_distance_tb <- 
  results_gr[queryHits(nearest_fli)] %>% 
  as_tibble(.) %>% 
  dplyr::mutate(nearest_insertion_distance = mcols(nearest_fli)$distance) %>% 
  tidyr::unite(coordinates.Siomi, seqnames, start, end, sep = " ") %>% 
  dplyr::select(-width) %>% 
  dplyr::mutate(gene_description = str_remove(gene_description, " \\[.*"))

# create FLIs table
flis_distance_tb <- 
  fli_gr[subjectHits(nearest_fli)] %>% 
  as_tibble(.) %>% 
  tidyr::unite(col = nearest_insertion_coord, seqnames, start, end, strand, sep = " ") %>% 
  dplyr::mutate(repName = repName %>% str_split(., "/") %>% purrr::map(., unique) %>% purrr::map(., str_c, collapse = "/") %>% unlist(.)) %>% 
  dplyr::select(nearest_insertion_coord, repName, repClass, insertion_class) 

# join together
results_insertions_tb <- cbind(results_distance_tb, flis_distance_tb) 

# save
readr::write_csv(results_insertions_tb, path = file.path(outpath, str_c("oocyte_Mov10l1.KO_vs_WT.upregulated_genes.distance_to_nearest_L1_IAP.csv")))