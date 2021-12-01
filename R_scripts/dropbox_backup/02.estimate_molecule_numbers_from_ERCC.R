### INFO: 
### DATE: Fri Nov 29 12:17:12 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Wu_2019_unpub_GSE133748/Analysis/ERCC_spike_expression")

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

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# calculate standard error of the mean FPKM value
standard.error <- function(x) {
  sqrt(var(x) / length(x))
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 93

# fpkm path
fpkm_path <- file.path(inpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.ERCC.FPKM.csv")

# ERCC info path
ercc_path <- "/common/DB/genome_reference/spike_in/ERCC92/ERCC92.info.txt"

######################################################## READ DATA
# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

# read ERCC info
ercc_tb <- readr::read_delim(ercc_path, delim = "\t")

######################################################## MAIN CODE
# clean ERCC info and calculate expected number of molecules
ercc_molecules <- 
  ercc_tb %>% 
  set_colnames(., c("num", "id", "subgroup", "conc_mix1", "conc_mix2", "expected_fc", "log2_mix1_mix2")) %>% # concentraction is in attomoles/ul
  dplyr::select(id, conc_mix1) %>% 
  dplyr::mutate(ercc_molecules_n = 
                  (conc_mix1 / 100000) * # mix was diluted 10^5 fold
                  1 *                    # 1 µl of the diluted ERCC mix was added to each sample
                  1/10^18 *              # Number of attomoles in a mole
                  6.02214179e23)         # Number of molecules in a mole

# join with ERCC FPKM
ercc_fpkm <- 
  ercc_molecules %>% 
  dplyr::select(-conc_mix1) %>%
  dplyr::left_join(., fpkm_tb %>% dplyr::select_at(vars(gene_id, starts_with("s_"))), by = c("id" = "gene_id")) %>% 
  tidyr::pivot_longer(-c(id, ercc_molecules_n),
                      names_to = "sample_id", 
                      values_to = "fpkm") %>% 
  dplyr::filter(ercc_molecules_n >= 1, 
                fpkm > 0) %>% 
  dplyr::select(gene_name = id, sample_id, fpkm, molecules_n = ercc_molecules_n)

# build linear regression model
linearMod <- lm(molecules_n ~ fpkm, data = ercc_fpkm)  # build linear regression model on full data
summary(linearMod)

### get FPKM of chosen genes, calculate number of molecules, join with ERCC, plot together 
# filter FPKM dataset, calculate number of molecules based on the linear regression model of ERCC spike-in
fpkm_tidy <- 
  fpkm_tb %>% 
  dplyr::select(gene_id, starts_with("s_")) %>% 
  tidyr::pivot_longer(-c(gene_id),
                      names_to = "sample_id", 
                      values_to = "fpkm") %>% 
  dplyr::mutate(molecules_n = round(linearMod$coefficients[1] + (linearMod$coefficients[2] * fpkm)))

## get and save wide tables
# number of molecules
fpkm_tidy %>%
  pivot_wider(-fpkm,
              names_from = "sample_id",
              values_from = "molecules_n") %T>%
  readr::write_csv(., file.path(outpath, "ensembl.93.GRCm38.p6.20180919.UCSCseqnames.ERCC.Wu_2019.n_molecules.csv"))

