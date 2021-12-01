### INFO: get expression of all genes in Fugaku's data
### DATE: Thu Oct 25 11:13:50 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(SummarizedExperiment)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### directory paths
# set experiment path
experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/ENCODE_2014_Nature_GSE49417"

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# set expression analysis path
expression_path <- file.path(experiment_path, "Analysis/expression")

# set documentation path
documentation_path <- file.path(experiment_path, "Data/Documentation")

# set working dir to expression analyis path 
setwd(expression_path)

# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()


### files paths
# summarizedOverlaps path
se_path <- list.files(expression_path, "ensembl.93.*.se.RDS", full.names = T)

# sample table path
sample_table_path <- list.files(documentation_path, ".*sample_table.csv", full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_dir, pattern = "ensembl.93.*UCSCseqnames.reducedExons.RDS$", full.names = T)

######################################################## READ DATA
# read summarizedExperiment from RDS file
se <- readRDS(file = se_path)

# read sample table
sample_table <- data.table::fread(file = sample_table_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  data.table(gene_id = names(.), width = .)

# set tissue order
tissue_order <- c("cns.E11.5", "cns.E14", "cns.E18", "frontallobe", "cortex", "cerebellum",
                  "stomach", "liver", "duodenum", "smintestine", "lgintestine", "colon", 
                  "lung", "heart", "bladder", "kidney", "thymus", "mammarygland", "spleen", 
                  "ovary", "testis", "placenta")

# get FPKM
fpkm_dt <-
  assay(se) %>%
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var = "gene_id") %>%
  set_colnames(., str_replace(colnames(.), ".Aligned.sortedByCoord.out.bam|.total.bam", "")) %>%
  data.table::as.data.table(.) %>% 
  data.table::melt(., 
                   id.vars = c("gene_id"),
                   variable.name = "sample_id", 
                   value.name = "counts") %>% 
  .[sample_table[, -"bam_path"], on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)), 
                                                        sample_id = i.tissue)] %>% 
  .[exons_width, on = "gene_id", fpkm := (fpm / round(width / 1E3, 3))] %>%
  .[, list(fpkm = round(mean(fpkm), 3)), by = .(gene_id, sample_id)] %T>% 
  # dcast(., gene_id ~ tissue, value.var = "fpkm") %>% 
  # setcolorder(., c("gene_id", tissue_order)) %T>% 
  readr::write_csv(., path = file.path(outpath, basename(se_path) %>% str_replace(., ".se.RDS", ".FPKM.csv")))
  
