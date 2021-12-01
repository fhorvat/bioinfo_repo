### INFO: 
### DATE: Fri Nov 02 19:59:36 2018
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

library(Biobase)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# set experiment path
exp_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Su_2004_ProcNatlAcadSciUSA_GSE1133"

# sample table path
sample_table_path <- list.files(exp_path, "*.sampleTable.csv", full.names = T)

# GCRMA normalized .RDS path
gcrma_path <- file.path(exp_path, "accessory_files", "GPL1073.gcrma_affy.RDS")

# mas5 normalized .RDS path
mas5_path <- file.path(exp_path, "accessory_files", "GPL1073.mas5_affy.RDS")

# probe info path
probe_info_path <- file.path(exp_path, "accessory_files/gnGNF1Musa_Mm_ENSG_23.0.0", "gnGNF1Musa_Mm_ENSG_mapping.txt")

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.geneInfo.csv$", full.names = T)

# accessory data path
accessory_data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data"

# hsapiens homology info
homology_info_path <- file.path(accessory_data_path, "ensembl.93.mouse_human.one2one_homologs.csv")

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# # read GCRMA .RDS
# gcrma_affy <- readRDS(gcrma_path)

# read Mas5 .RDS
mas5_affy <- readRDS(mas5_path)

# read probe info
probe_info <- readr::read_delim(probe_info_path, delim = "\t") 

# read info about chosen genes
genes_info <- readr::read_csv(file = genes_info_path)

# # read info about one2one homologs between mouse and human
homologs_info <- readr::read_csv(file = homology_info_path)

######################################################## MAIN CODE
# get gene_id of protein coding genes which have only one2one homology in humans
filtered_genes <- 
  homologs_info %>% 
  dplyr::select(gene_id) %>% 
  left_join(., genes_info, by = "gene_id") %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# clean probe info
probe_info %<>% 
  dplyr::select(gene_id = `Probe Set Name`, probe_id = `Affy Probe Set Name`) %>% 
  dplyr::distinct(.) %>% 
  dplyr::mutate(gene_id = str_remove(gene_id, "_at")) %>% 
  dplyr::filter(str_detect(probe_id, "_.{1}_at$")) %>% 
  # dplyr::filter(gene_id %in% filtered_genes) %>%
  as.data.table(.)

## get expression data
# set tissue
tissue <- "skeletalmuscle"

# get binned expression in tissue
binned_exp <- 
  Biobase::exprs(mas5_affy) %>% 
  as.data.table(., keep.rownames = "probe_id") %>% 
  .[probe_info, on = "probe_id", gene_id := i.gene_id] %>% 
  .[!(is.na(gene_id)), -"probe_id"] %>%
  melt(., id.vars = "gene_id", variable.name = "geo_accession", value.name = "exp_value", variable.factor = F) %>% 
  .[, .(exp_value = round(mean(exp_value), 3)), by = list(gene_id, geo_accession)] %>%
  .[sample_table, on = "geo_accession", sample_tissue := i.sample_tissue] %>% 
  .[, .(exp_value = round(mean(exp_value), 3)), by = list(gene_id, sample_tissue)] %>% 
  .[order(gene_id, exp_value)] %>% 
  .[, position := order(exp_value), by = "gene_id"] %>% 
  .[sample_tissue == tissue & exp_value > 0] %>% 
  .[order(position), bin_relative := as.numeric(Hmisc::cut2(1:.N, g = 50))] %>% 
  .[order(exp_value), bin_absolute := as.numeric(Hmisc::cut2(exp_value, g = 50))] %>% 
  .[order(exp_value)]

# save binned FPKMs
# readr::write_csv(binned_exp, file.path(outpath, str_c("grid.Su.mas5.", tissue, ".50.bins.protein_coding.human_homologs.csv")))
readr::write_csv(binned_exp, file.path(outpath, str_c("grid.Su.mas5.", tissue, ".50.bins.csv")))

  