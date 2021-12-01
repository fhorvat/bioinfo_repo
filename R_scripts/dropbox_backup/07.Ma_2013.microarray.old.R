### INFO: read Ma 2013 microarray data
### DATE: Mon Mar 12 23:39:34 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Ma_2013_BiolReprod_GSE27049")

######################################################## LIBRARIES
library(affy)
library(gcrma)

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

# cel files path
cel_path <- list.files(inpath, pattern = ".CEL.gz", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- read.csv(file = file.path(inpath, "Ma_2013_BiolReprod_GSE27049.sample_table.csv"), stringsAsFactors = F)

######################################################## MAIN CODE
# set rownames
rownames(sample_table) <- sample_table$sample_name
sample_table$sample_name <- NULL

# set sample names
sample_names <- sample_table$array_id
names(sample_names) <- rownames(sample_table)

# read affymetrix microarray set
rawAffy <- affy::read.affybatch(filenames = sample_table$cel_path, phenoData = Biobase::AnnotatedDataFrame(sample_table[, 1, drop = F]))

# normalize using GCRMA algorithm
gcrma_affy <- gcrma::gcrma(rawAffy)

# extract expression values, remove those for which one or more controls has intesity < 10 
exprs_matrix <- Biobase::exprs(gcrma_affy)
exprs_matrix <- 2^exprs_matrix

# test with SAM algorithm which probes are significantly up- or down-regulated
sam_out_sub <- samr::SAM(x = exprs_matrix, 
                         y = rep(c(1, 2), each = 3),
                         resp.type = "Two class unpaired",
                         nperms = 100, 
                         testStatistic = "standard", 
                         fdr.output	= 0.05, 
                         logged2 = F)

up_genes_sub <- sam_out_sub$siggenes.table$genes.up
