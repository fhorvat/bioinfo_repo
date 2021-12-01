### INFO: expression of lnc1 and other genes in ENCODE dataset
### DATE: Tue Aug 21 14:25:50 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path 
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = "ensembl.91.*[0-9]{6}.UCSCseqnames.reducedExons.RDS$", full.names = T)

# gene info path
gene_info_path <- list.files(genome_path, "ensembl.91.GRCm38.p5.*.UCSCseqnames.geneInfo.csv", full.names = T)

# sample table path
sample_table_path <- list.files(inpath, pattern = ".*.sample_table.csv", full.names = T)

# summarizedOverlaps path
se_path <- file.path(inpath, "ensembl.91.GRCm38.p5.20180512.lncKO.summarizedOveralaps.RDS")

######################################################## READ DATA
# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read info about all ENSEMBL annotated genes
gene_info <- readr::read_csv(gene_info_path)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read summarizedExperiment from RDS file
se <- readRDS(file = se_path)

######################################################## MAIN CODE
# get resequenced samples
reseq_samples <- 
  sample_table$sample_id %>% 
  .[str_detect(., "reseq")] %>% 
  str_remove(., ".reseq")

# remove resequenced samples from sample table
sample_table %<>%
  dplyr::filter(!(sample_id %in% reseq_samples)) %T>%
  readr::write_csv(., path = file.path(outpath, "lncKO.filtered.sample_table.csv"))

# get total length of all exons for each transcript
exons_width <-
  width(exons_gr) %>%
  sum(.) %>%
  tibble(gene_id = names(.), width = .)

# get data.frame of counts, transform to FPKM
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
  dplyr::right_join(., sample_table %>% dplyr::select(sample_id, genotype, library_size), by = "sample_id") %>%
  dplyr::left_join(., exons_width, by = "gene_id") %>%
  dplyr::mutate(library_size = round(library_size / 1E6, 6),
                width = round(width / 1E3, 3),
                fpm = (counts / library_size),
                fpkm = (fpm / width)) %>% 
  dplyr::select(gene_id, sample_id, genotype, fpkm)

# save individual samples FPKM
fpkm_df %>% 
  dplyr::select(-genotype) %>% 
  tidyr::spread(sample_id, fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, "lncRNA_KO.individual_samples.FPKM.csv"))

# save average genotype FPKM
fpkm_df %>% 
  dplyr::group_by(gene_id, genotype) %>% 
  dplyr::summarise(fpkm = mean(fpkm)) %>% 
  tidyr::spread(genotype, fpkm) %T>%
  readr::write_csv(., path = file.path(outpath, "lncRNA_KO.mean_genotype.FPKM.csv"))




