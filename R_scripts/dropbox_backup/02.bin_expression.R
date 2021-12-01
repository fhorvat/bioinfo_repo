### INFO: 
### DATE: Fri Nov 02 19:59:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/biogps")

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

# biogps GCRMA normalized table path
biogps_path <- file.path(exp_path, "accessory_files", "GNF1M_plus_macrophage_small.bioGPS.txt.gz")

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

# read biogps expression
biogps_exp <- data.table::fread(biogps_path)
data.table::setnames(biogps_exp, make.unique(colnames(biogps_exp)))
  
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

### bin expression in different tissue
# get average expression values in all tissues
biogps_mean <- 
  biogps_exp %>%
  data.table::setnames(., 1, "probe_id") %>%
  .[probe_info, on = "probe_id", gene_id := i.gene_id] %>% 
  .[!(is.na(gene_id)), -"probe_id"] %>%
  melt(., id.vars = "gene_id", variable.name = "sample_id", value.name = "exp_value", variable.factor = F) %>% 
  .[, .(exp_value = round(mean(exp_value), 3)), by = list(gene_id, sample_id)] %>%
  .[, sample_tissue := str_remove(sample_id, "\\.1$") %>% str_replace_all(., " ", "_")] %>% 
  .[, .(exp_value = round(mean(exp_value), 3)), by = list(gene_id, sample_tissue)] %>% 
  .[order(gene_id, exp_value)] %>% 
  .[, position := order(exp_value), by = "gene_id"] %>% 
  .[]

# set tissue list
tissue_list <- c("skeletalmuscle", "liver", "oocyte", "blastocysts", "embryonic_stem_no_feeder", "embryonic_stem_feeder_layer", "fertilizedegg")

# set number of bins
bin_number <- 50

# loop through tissues
for(tissue in tissue_list){
  
  binned_exp <- 
    biogps_mean %>% 
    .[sample_tissue == tissue & exp_value > 0] %>% 
    .[order(position), bin_relative := as.numeric(Hmisc::cut2(1:.N, g = bin_number))] %>% 
    .[order(exp_value), bin_absolute := as.numeric(Hmisc::cut2(exp_value, g = bin_number))] %>% 
    .[order(exp_value)]
  
  # save as .RDS
  saveRDS(binned_exp, file.path(outpath, str_c("grid.biogps.gcrma", tissue, bin_number, "bins.RDS", sep = ".")))
  
}

