### INFO: 
### DATE: Thu Oct 25 10:46:52 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Su_2004_ProcNatlAcadSciUSA_GSE1133")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(AnnotationForge)
library(makecdfenv)
library(affy)
library(gcrma)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# cel files
cel_path <- list.files(inpath, "*.CEL.gz", full.names = T)

# raw array expression path
exp_mat_path <- file.path(inpath, "accessory_files", "GSE1133-GPL1073_series_matrix.txt.gz")

######################################################## READ DATA
# raw_exp <- readr::read_delim(file = exp_mat_path, skip = 70, col_names = T, delim = "\t")

# read info about array samples
sample_table <- 
  readr::read_delim(file = exp_mat_path, skip = 44, n_max = 2, col_names = F, delim = "\t") %>% 
  tidyr::gather(key = "variable", value = "geo", -X1) %>% 
  tidyr::spread(key = "X1", value = "geo") %>% 
  dplyr::select(sample_id = `!Sample_title`, geo_accession = `!Sample_geo_accession`) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "M[A-Z]{3}\\d{9}"), 
                replicate = str_extract(sample_id, "^[A,B]{1}"), 
                replicate = ifelse(replicate == "A", "r1", "r2"), 
                sample_id = str_remove(sample_id, "^[A,B]{1}"), 
                sample_tissue = sample_id,
                sample_id = str_c(sample_id, "_", replicate)) %>% 
  dplyr::arrange(sample_id) %>% 
  dplyr::left_join(., tibble(geo_accession = (basename(cel_path) %>% str_remove(., ".CEL.gz")), 
                             raw_path = cel_path), 
                   by = "geo_accession") %T>% 
  readr::write_csv(., path = file.path(outpath, "Su_2004_ProcNatlAcadSciUSA_GSE1133.sampleTable.csv"))

######################################################## MAIN CODE
# ##### 
# ## CDF
# # build cdf package package, manually change DESCRIPTION file to proper format and install package "using R CMD INSTALL gngnf1musacdf"
# makecdfenv::make.cdf.package(filename = "GPL1073.CDF",
#                              packagename  = "gngnf1musacdf",
#                              cdf.path = outpath,
#                              package.path = outpath,
#                              compress = FALSE,
#                              author = "fhorvat",
#                              maintainer = "fhorvat",
#                              version = packageDescription("makecdfenv", fields ="Version"),
#                              species = "Mus_musculus",
#                              unlink = FALSE,
#                              verbose = F)
# 
# ## probe annotation
# # read and clean tab file
# readr::read_delim(file = file.path(outpath, "GPL1073.tab.gz"), delim = "\t") %>% 
#   dplyr::select(-`Serial Order`) %T>% 
#   readr::write_delim(., path = file.path(outpath, "GPL1073.clean.tab"), col_names = T, delim = "\t")
# 
# # build probe annotation package
# # manually change DESCRIPTION file to proper format and install package using "R CMD INSTALL gngnf1musaprobe"
# AnnotationForge::makeProbePackage(arraytype = "gngnf1musaprobe", 
#                                   datafile = file.path(inpath, "GPL1073.clean.tab"),
#                                   importfun = "getProbeDataAffy", 
#                                   maintainer = c("fhorvat", "fihorvat@gmail.com"), 
#                                   version = "0.0.1", 
#                                   species = "Mus musculus", 
#                                   pkgname = "gngnf1musaprobe", 
#                                   outdir  = ".", 
#                                   quiet = FALSE, 
#                                   check = TRUE, 
#                                   build = TRUE, 
#                                   unlink = TRUE)

#####
## read and normalize microarrays
sample_table %<>%
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "geo_accession")

# read raw files
rawAffy <- affy::read.affybatch(filenames = sample_table$raw_path, 
                                phenoData = AnnotatedDataFrame(sample_table))

# raw expression
raw_exp <- 
  Biobase::exprs(rawAffy) %>% 
  as.data.table(., keep.rownames = "probe_id")

# mas5 normalize
mas5_affy <- affy::mas5(rawAffy)

# gcrma normalize
gcrma_affy <- gcrma::gcrma(rawAffy)

# save as .RDS
saveRDS(gcrma_affy, "GPL1073.gcrma_affy.RDS")
saveRDS(mas5_affy, "GPL1073.mas5_affy.RDS")



