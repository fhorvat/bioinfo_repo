### INFO: creates table with relations between ensembl and USCS seqnames using assembly report from NCBI
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/DB/genome_reference/golden_hamster/MesAur1.0.GCA_000349665.1")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(biomaRt)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath
inpath <- getwd()

# set outpath
outpath <- getwd()

# ensembl gtf path
ensembl_gtf_path <- file.path(inpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.gtf.gz")

# ncbi gtf path
ncbi_gtf_path <- file.path(inpath, "GCF_000349665.1_MesAur1.0_genomic.gtf")

# changed last column path
last_column_path <- file.path(inpath, "GCF_000349665.1_MesAur1.0_genomic.PIWIL3.last_column_ensembl_style.txt")

# gene info path
gene_info_path <- file.path(inpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.geneInfo.csv")

######################################################## READ DATA
# read gtfs
ensembl_gtf <- read_delim(file = ensembl_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))
ncbi_gtf <- read_delim(file = ncbi_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

# read last column
last_column <- readr::read_lines(last_column_path)

# read gene info
gene_info <- readr::read_csv(gene_info_path)

######################################################## MAIN CODE
### add PIWIL3 from NCBI to .gtf
# find gene in NCBI gtf, add last column in ensembl version
ncbi_gtf_filt <- 
  ncbi_gtf %>% 
  dplyr::filter(str_detect(X9, "rna20892|gene14314")) %>% 
  dplyr::mutate(X1 = "KB708214.1", 
                X3 = replace(X3, X3 == "mRNA", "transcript"), 
                X9 = last_column)

# bind with original gtf
ensembl_gtf_piwil3 <- bind_rows(ensembl_gtf, ncbi_gtf_filt)

# save
write.table(x = ensembl_gtf_piwil3, file = file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.NCBI_PIWIL3.gtf"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


### add info about PIWIL3 to gene info table
# create table
piwil3_info <- tibble(gene_id = "gene14314", 
                      seqnames = "KB708214.1", 
                      start = 1860221, 
                      end = 1889111, 
                      strand = "+", 
                      gene_name = "Piwil3", 
                      gene_biotype = "protein_coding",
                      gene_description = "piwi like RNA-mediated gene silencing 3 [Source:NCBI gene;Acc:101844316]")

# bind with gene info
gene_info_piwil3 <- bind_rows(gene_info, piwil3_info)
  
# save
readr::write_csv(gene_info_piwil3, file.path(outpath, "ensembl.99.MesAur1.0.20200415.UCSCseqnames.NCBI_PIWIL3.geneInfo.csv"))
  
  