### INFO: read Ivanova 2017 YTHDF2 KO microarray data
### DATE: Mon Apr 09 23:25:24 2018
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
library(mogene20stprobeset.db)
library(mogene20sttranscriptcluster.db)
library(genefilter)
library(limma)
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
rma_affy <- oligo::rma(rawAffy, background = T, normalize = T)

# MII only
rma_affy <- rma_affy[, str_detect(sampleNames(rma_affy), "mii")]

# limma
# MII only
data.matrix <- exprs(rma_affy)
ph <- rma_affy@phenoData
groups <- ph@data$genotype
f <- factor(groups, levels = c("KO", "WT"))
design <- model.matrix(~ 0 + f)
colnames(design) <- c("KO", "WT")
data.fit <- lmFit(data.matrix, design)
contrast.matrix <- makeContrasts(KO-WT, levels = design)
data.fit.con <- contrasts.fit(data.fit, contrast.matrix)
data.fit.eb <- eBayes(data.fit.con)
top <- topTable(data.fit.eb, n = 100000, adjust.method = "BH", coef = "KO - WT")
top_filt <- top[top$adj.P.Val <= 0.05, ]
top_filt <- top_filt[abs(top_filt$logFC) >= 2, ]

# get upregulated gens with logFC > 2, save
up_genes <-
  top_filt %>%
  tibble::rownames_to_column(., var = "probe_id") %>% 
  as.tibble(.) %>%
  dplyr::mutate(logFC = as.numeric(logFC)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  dplyr::filter(logFC >= 2) %>%
  dplyr::mutate(gene_id1 = mapIds(mogene20sttranscriptcluster.db, keys = probe_id, column = "ENSEMBL", keytype = "PROBEID", multiVals = "first"), 
                gene_id2 = mapIds(mogene20stprobeset.db, keys = probe_id, column = "ENSEMBL", keytype = "PROBEID", multiVals = "first")) %>% 
  mutate_cond(!is.na(gene_id1), gene_id = gene_id1) %>% 
  mutate_cond(!is.na(gene_id2), gene_id = gene_id2) %>% 
  dplyr::select(gene_id, logFC) %T>% 
  readr::write_csv(., path = file.path(outpath, str_c("affy.MoGene2.0st.", experiment_name, ".KOvsWT.upregulated.csv")))

