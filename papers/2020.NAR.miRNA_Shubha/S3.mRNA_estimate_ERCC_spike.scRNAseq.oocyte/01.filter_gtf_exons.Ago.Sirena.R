### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/ERCC_estimates")

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
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.geneInfo\\.csv$"), full.names = T)

# gtf path
gtf_path <- list.files(path = genome_dir, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames\\.gtf\\.gz$"), full.names = T)

######################################################## READ DATA
# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read gtf
gtf_tb <- readr::read_delim(file = gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c"))

######################################################## MAIN CODE
### get only exons and genes from .gtf
gtf_filt <- 
  gtf_tb %>% 
  dplyr::filter(X3 %in% c("exon", "gene"))


### get Ago exons
# get gene ID's of Ago genes
gene_ids_list <- 
  genes_info %>% 
  dplyr::filter(str_detect(gene_name, "Ago[1-4]")) %$%
  gene_id

# get last 8 exons for each Ago protein 
Ago_filt <- 
  gtf_filt %>% 
  dplyr::filter(str_detect(X9, str_c(gene_ids_list, collapse = "|"))) %>% 
  dplyr::mutate(gene_id = str_extract(X9, str_c(gene_ids_list, collapse = "|")), 
                transcript_index = str_extract(X9, '(?<=transcript_name ")Ago[1-4]-[0-9]{3}') %>% 
                  str_remove(., "Ago[1-4]-") %>% 
                  as.numeric(.), 
                exon_number = str_extract(X9, '(?<=exon_number ")[1-9]+') %>% 
                  as.numeric(.)) %>% 
  dplyr::filter(transcript_index == 201 | X3 == "gene") %>% 
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(total_exons = max(na.omit(exon_number))) %>% 
  dplyr::filter((exon_number >= total_exons - 7) | (X3 == "gene")) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(X1:X9)


### get Sirena1 exons
Sirena1_filt <- 
  gtf_filt %>% 
  dplyr::filter(str_detect(X9, "C86187"),
                X3 %in% c("exon", "gene")) %>% 
  dplyr::filter(str_detect(X9, "C86187-202") | X3 == "gene") %>% 
  dplyr::filter(str_detect(X9, 'exon_number \"1\"|exon_number \"4\"') | X3 == "gene")


### remove normal Ago/Sirena, add new Ago/Sirena
gtf_filt_fixed <- 
  gtf_filt %>% 
  dplyr::filter(!str_detect(X9, str_c(c(gene_ids_list, "C86187"), collapse = "|"))) %>% 
  rbind(., Ago_filt, Sirena1_filt) %>% 
  arrange(X1, X9)

# save as .gtf
gtf_name <- basename(gtf_path) %>% str_replace(., "\\.gtf\\.gz", ".Ago.Sirena.short.gtf")
write.table(x = gtf_filt_fixed, file = file.path(outpath, gtf_name), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


 