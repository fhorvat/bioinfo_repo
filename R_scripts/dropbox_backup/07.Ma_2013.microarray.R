### INFO: read Ma 2013 microarray data
### DATE: Mon Mar 12 23:39:34 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Ma_2013_BiolReprod_GSE27049")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(affy)
library(limma)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# sample table path
sample_table_path <- file.path(inpath, "Ma_2013_BiolReprod_GSE27049.sample_table.csv")

# cel files path
cel_path <- list.files(inpath, pattern = ".CEL.gz", full.names = F)

######################################################## READ DATA
# read sample table
sample_table <- read_csv(sample_table_path)

######################################################## MAIN CODE
# add path to cell table, prepare for AnnotatedDataFrame
sample_table %<>% 
  dplyr::left_join(., tibble(cel_path, sample_id = str_replace(cel_path, ".CEL.gz", "")), by = "sample_id") %>% 
  as.data.frame(.) %>% 
  column_to_rownames(var = "sample_id")

# import "phenotype" data describing the experimental design
phenoData <- AnnotatedDataFrame(sample_table[, 1, drop = F])

# RMA normalization
eset <- affy::justRMA(filenames = sample_table$cel_path, phenoData = phenoData, compress = T)

# differential expression
combn <- factor(pData(phenoData)[,1])
design <- model.matrix(~combn)

# fit each probeset to model
fit <- lmFit(eset, design)  

# empirical Bayes adjustment
efit <- eBayes(fit)

# table of differentially expressed probesets
results_df <- 
  topTable(efit, coef = 2, n = Inf) %>% 
  tibble::rownames_to_column(., var = "affy_probe") %>% 
  as.tibble(.) %>% 
  dplyr::filter(adj.P.Val < 0.05, 
                logFC > 1.5) %>%
  dplyr::arrange(desc(logFC))

# Extract the table from experiment
out <- extract_tables("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4434936/bin/supp_112.105312_105312SupData.pdf")

# bind to data.frame
pdf_df <- 
  lapply(out, as.data.frame) %>% 
  dplyr::rbind_all(.)