### INFO: read Ma 2013 microarray data
### DATE: Mon Mar 12 23:39:34 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Ivanova_2017_MolCell_EMTAB5576")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(oligo)
library(gcrma)
library(Biobase)
library(pd.mogene.2.0.st)
library(genefilter)
library(samr)
library(biomaRt)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set experiment name
experiment_name <- "Ivanova2017"

# set main outpath
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis"

# sample table path
sample_table_path <- file.path(inpath, "Ivanova_2017_MolCell_EMTAB5576.sampleTable.txt")

# cel files path
cel_path <- list.files(inpath, pattern = ".CEL", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- read_delim(sample_table_path, delim = "\t")

######################################################## MAIN CODE
# add path to cell table, prepare for AnnotatedDataFrame
sample_table_filt <- 
  sample_table %>% 
  dplyr::select(sample_name = `Source Name`, stage = `Characteristics[cell type]`, array_id = `Array Data File`) %>% 
  dplyr::mutate(genotype = toupper(stringr::str_extract(sample_name, "ko|ctrl")), 
                genotype = replace(genotype, genotype == "CTRL", "WT"), 
                cel_path = file.path(inpath, array_id), 
                stage = stringr::str_replace(stage, " Oocyte", "")) %>% 
  as.data.frame(.) %>% 
  column_to_rownames(var = "sample_name")

# read affymetrix microarray set
rawAffy <- oligo::read.celfiles(filenames = sample_table_filt$cel_path, 
                                phenoData = AnnotatedDataFrame(sample_table_filt %>% dplyr::select(stage, genotype)))

# normalize using RMA algorithm
rma_affy <- oligo::rma(rawAffy)

# extract expression values
exprs_matrix <-
  Biobase::exprs(rma_affy) %>%
  .[, str_detect(colnames(.), "mii")] %>%
  2^.

# test with SAM algorithm which probes are significantly up- or down-regulated (control = 1, KO = 2)
sam_out_sub <- samr::SAM(x = exprs_matrix, 
                         y = rep(c(1, 2), times = c(4, 3)),
                         geneid = rownames(exprs_matrix), 
                         resp.type = "Two class unpaired",
                         nperms = 100, 
                         testStatistic = "standard", 
                         fdr.output	= 0.05, 
                         logged2 = F)

# get upregulated gens with logFC > 2, save
up_genes <-
  sam_out_sub$siggenes.table$genes.up %>%
  as.tibble(.) %>%
  dplyr::select(probe_id = `Gene Name`, logFC = `Fold Change`) %>% 
  dplyr::mutate(logFC = as.numeric(logFC)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  dplyr::filter(abs(logFC) >= 2) %>% 
  dplyr::mutate(gene_id1 = mapIds(mogene20sttranscriptcluster.db, keys = probe_id, column = "ENSEMBL", keytype = "PROBEID", multiVals = "first"), 
                gene_id2 = mapIds(mogene20stprobeset.db, keys = probe_id, column = "ENSEMBL", keytype = "PROBEID", multiVals = "first")) %>% 
  mutate_cond(!is.na(gene_id1), gene_id = gene_id1) %>% 
  mutate_cond(!is.na(gene_id2), gene_id = gene_id2) %>% 
  dplyr::select(gene_id, logFC) %T>% 
  readr::write_csv(., path = file.path(outpath, str_c("affy.MoGene2.0st.", experiment_name, ".KOvsWT.upregulated.csv")))

