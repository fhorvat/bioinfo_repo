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
### main paths
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


### array paths
# set experiment path
array_exp_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Su_2004_ProcNatlAcadSciUSA_GSE1133"

# sample table path
sample_table_path <- list.files(array_exp_path, "*.sampleTable.csv", full.names = T)

# mas5 normalized .RDS path
mas5_path <- file.path(array_exp_path, "accessory_files", "GPL1073.mas5_affy.RDS")

# probe info path
probe_info_path <- file.path(array_exp_path, "accessory_files/gnGNF1Musa_Mm_ENSG_23.0.0", "gnGNF1Musa_Mm_ENSG_mapping.txt")


### ngs paths
# ENCODE FPKM path
encode_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.FPKM.csv"

# Fugaku FPKM path
fugaku_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Fugaku_RNAseq/Analysis/expression/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.Fugaku.FPKM.csv"


### gene info paths
# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.geneInfo.csv$", full.names = T)

# accessory data path
accessory_data_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/expression_grids/accessory_data"

# hsapiens homology info
homology_info_path <- file.path(accessory_data_path, "ensembl.93.mouse_human.one2one_homologs.csv")

######################################################## READ DATA
### array data
# read sample table
sample_table <- data.table::fread(sample_table_path)

# read Mas5 .RDS
mas5_affy <- readRDS(mas5_path)

# read probe info
probe_info <- readr::read_delim(probe_info_path, delim = "\t") 


### NGS data
# read ENCODE FPKM table
encode_fpkm <- data.table::fread(file = encode_path)

# read Fugaku FPKM table
fugaku_fpkm <- data.table::fread(file = fugaku_path)


### gene info data
# read info about chosen genes
genes_info <- readr::read_csv(file = genes_info_path)

# read info about one2one homologs between mouse and human
homologs_info <- readr::read_csv(file = homology_info_path)

######################################################## MAIN CODE
#### CLEAN DATA ####
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
  dplyr::filter(gene_id %in% filtered_genes, 
                str_detect(probe_id, "_.{1}_at$")) %>% 
  as.data.table(.)

#### ARRAY ####
# set tissue
tissue <- "oocyte"

# get expressed genes on array in oocyte
array_exp <- 
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
  .[sample_tissue == tissue & exp_value > 0]

# get expressed genes on NGS in oocyte
fugaku_fpkm_clean <- fugaku_fpkm[(gene_id %in% array_exp$gene_id) & (sample_id == "s_GV.WE") & (fpkm > 0), ]

# filter out 
array_exp_binned <- 
  array_exp %>% 
  .[gene_id %in% fugaku_fpkm_clean$gene_id, ] %>%
  .[order(position), bin_relative := as.numeric(Hmisc::cut2(1:.N, g = 50))] %>% 
  .[order(exp_value), bin_absolute := as.numeric(Hmisc::cut2(exp_value, g = 50))] %>% 
  .[order(exp_value)]

# save binned FPKMs
readr::write_csv(array_exp_binned, file.path(outpath, str_c("grid.Su.mas5.", tissue, ".50.bins.filtered.csv")))


#### NGS ####
##### ENCODE ####
# get FPKM's for only protein coding genes, get mean expression across all tissue
encode_fpkm_clean <- encode_fpkm[sample_id != "ovary", ]

# get FPKM's for only protein coding genes in oocyte, bin FPKM to 61 bins of equal size
fugaku_fpkm_relative <- 
  fugaku_fpkm[sample_id == "s_GV.WE", ] %>% 
  .[, sample_id := "oocyte"] %>%
  rbind(., encode_fpkm_clean) %>% 
  .[order(gene_id, fpkm)] %>% 
  .[, position := order(fpkm), by = "gene_id"] %>% 
  .[sample_id == "oocyte", ] %>%
  .[order(position), bin_relative := as.numeric(Hmisc::cut2(1:.N, g = 50))] %>% 
  .[order(position)]

### Fugaku ####
# get FPKM's for only protein coding genes in oocyte, bin FPKM to 61 bins of equal size
fugaku_fpkm_absolute <- 
  fugaku_fpkm[sample_id == "s_GV.WE", ] %>% 
  .[, sample_id := "oocyte"] %>%
  .[gene_id %in% array_exp_binned$gene_id, ] %>% 
  .[order(fpkm), bin_absolute := as.numeric(Hmisc::cut2(fpkm, g = 50))] %>% 
  .[fugaku_fpkm_relative, on = "gene_id", `:=`(bin_relative = i.bin_relative)] %>% 
  .[]

# save binned FPKMs
readr::write_csv(fugaku_fpkm_absolute, file.path(outpath, "grid.Fugaku.Encode.oocyte.50.bins.filtered.csv"))


