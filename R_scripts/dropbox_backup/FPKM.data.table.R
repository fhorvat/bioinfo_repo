### INFO: 
### DATE: Thu Oct 25 21:44:54 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
# options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/other_projects/lab_meetings_Zagreb/2018_10_16")

######################################################## LIBRARIES
library(SummarizedExperiment)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)
library(tibble)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
se <- readRDS(file = "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.test.RDS")
colnames(se) <- str_remove(colnames(se), "\\.Aligned.sortedByCoord.out.bam")

sample_table <- data.table::fread("sample_table.csv")
exons_table <- data.table::fread("exons_width.csv")

######################################################## MAIN CODE
### get FPKM from counts
fpkm_df <-
  assay(se) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  data.table::as.data.table(.)

# melt to long format
fpkm_df <-  data.table::melt(fpkm_df, id.vars = c("gene_id"), variable.name = "sample_id", value.name = "counts")

# join with sample table
fpkm_df[sample_table[, -"bam_path"], on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)),
                                                            genotype = i.genotype,
                                                            stage = i.stage)]

# join with exons table
fpkm_df[exons_table, on = "gene_id", fpkm := (fpm / round(width / 1E3, 3))]

# get mean FPKM across developmental stage / genotype combination
fpkm_df <- fpkm_df[, list(avg_fpkm = round(mean(fpkm), 3)), by = .(gene_id, stage, genotype)]

# cast to wide format
fpkm_df <- dcast(fpkm_df, gene_id ~ stage + genotype, value.var = "avg_fpkm")

### use magrittr pipes 
fpkm_df <-
  assay(se) %>%
  as.data.table(., keep.rownames = "gene_id") %>% 
  melt(., 
       id.vars = c("gene_id"),
       variable.name = "sample_id", 
       value.name = "counts") %>% 
  .[sample_table[, -"bam_path"], on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)), 
                                                        genotype = i.genotype, 
                                                        stage = i.stage)] %>% 
  .[exons_table, on = "gene_id", fpkm := (fpm / round(width / 1E3, 3))] %>%
  .[, list(avg_fpkm = round(mean(fpkm), 3)), by = .(gene_id, stage, genotype)] %>% 
  dcast(., gene_id ~ stage + genotype, value.var = "avg_fpkm")

# use data.table pipes
fpkm_df <- 
  as.data.table(assay(se), keep.rownames = "gene_id")[
    , melt(.SD, id.vars = c("gene_id"), variable.name = "sample_id", value.name = "counts") ][
      sample_table[, -"bam_path"], on = "sample_id", `:=`(fpm = (counts / round(library_size / 1E6, 6)), genotype = i.genotype, stage = i.stage) ][
        exons_table, on = "gene_id", fpkm := (fpm / round(width / 1E3, 3)) ][
          , .(avg_fpkm = round(mean(fpkm), 3)), by = .(gene_id, stage, genotype) ][
            , dcast(.SD, gene_id ~ stage + genotype, value.var = "avg_fpkm")]
