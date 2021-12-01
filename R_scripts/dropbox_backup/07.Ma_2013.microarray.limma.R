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
library(gcrma)
library(Biobase)
library(mouse4302.db)
library(genefilter)
library(samr)
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
cel_path <- list.files(inpath, pattern = ".CEL.gz", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- read_csv(sample_table_path)

######################################################## MAIN CODE
# add path to cell table, prepare for AnnotatedDataFrame
sample_table %<>% 
  as.data.frame(.) %>% 
  column_to_rownames(var = "sample_name")

# set sample names
sample_names <- 
  sample_table$array_id %>% 
  magrittr::set_names(rownames(sample_table))

# read affymetrix microarray set
rawAffy <- affy::read.affybatch(filenames = sample_table$cel_path,
                                phenoData = AnnotatedDataFrame(sample_table[, 1, drop = F]))

# normalize using RMA algorithm
rma_affy <- rma(rawAffy)

# create design matrix
f <- factor(sample_table$treatment, levels = c("control", "Dcp1a", "Dcp2", "Dcp1a_Dcp2"))
design <- model.matrix(~0 + f)
colnames(design) <- c("control", "Dcp1a", "Dcp2", "Dcp1a_Dcp2")

# To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(rma_affy, design)
contrast.matrix <- makeContrasts(Dcp1a - control, Dcp2 - control, Dcp1a_Dcp2 - control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# compute the pooled sample variance for each gene and sample mean expression of each gene in each treatment
fitTrtMean <- lmFit(rma_affy, design_treat, cor = corfit$consensus.correlation)

# create the coefficient matrix for the contrasts
??duplicateCorrelation



# extract expression values, remove those for which one or more controls has intesity < 10 
exprs_matrix <- 
  Biobase::exprs(rma_affy) %>% 
  as.data.frame(.) %>% 
  magrittr::set_colnames(., value = sample_names[colnames(.)]) %>% 
  tibble::rownames_to_column(., var = "probe_id") %>% 
  as.tibble(.) %>% 
  dplyr::select_at(.vars = vars(probe_id, matches("control|Dcp1a_Dcp2"))) %>% 
  # dplyr::filter(rowSums(.[, 2:4] > log2(10)) >= 2) %>%
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "probe_id") %>% 
  as.matrix(.) %>% 
  2^.

# test with SAM algorithm which probes are significantly up- or down-regulated
sam_out_sub <- samr::SAM(x = exprs_matrix, 
                         y = rep(c(1, 2), each = 3),
                         geneid = rownames(exprs_matrix), 
                         resp.type = "Two class unpaired",
                         nperms = 100, 
                         testStatistic = "standard", 
                         fdr.output	= 0.05, 
                         logged2 = F)

# get upregulated gens with logFC > 1.5
up_genes <-
  sam_out_sub$siggenes.table$genes.up %>%
  as.tibble(.) %>%
  dplyr::select(probe_id = `Gene Name`, logFC = `Fold Change`) %>% 
  dplyr::mutate(logFC = as.numeric(logFC)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  dplyr::filter(logFC >= 1.5 ) %>% 
  dplyr::mutate(gene_id = mapIds(mouse4302.db, keys = probe_id, column = "ENSEMBL", keytype = "PROBEID", multiVals = "first"), 
                gene_symbol = mapIds(mouse4302.db, keys = probe_id, column = "SYMBOL", keytype = "PROBEID", multiVals = "first"))

