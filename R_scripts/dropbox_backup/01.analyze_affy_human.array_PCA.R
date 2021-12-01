### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/arrays/Joshi_2007_BMCDevBiol_GSE5558")

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

library(GEOquery)
library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

######################################################## READ DATA
### get data from GEO
# chip annotation
annotation <- getGEO("GPL3771", destdir = outpath, AnnotGPL = T)

# series matrix
eset <- getGEO("GSE5558", destdir = outpath)

######################################################## MAIN CODE
### series matrix from GEO includes only log2FC values and not individual array intensities,  
### so we need to get expression value for individual arrays
# get list of array expression values from GEO 
array_list <- 
  purrr::map(sampleNames(eset), getGEO, destdir = outpath) %>% 
  set_names(sampleNames(eset))

# get sample table
sample_tb <- purrr::map(names(array_list), function(array_name){
  
  # get one sample
  array_sample <- array_list[[array_name]] 
  
  # get RMA values
  sample_tb <- tibble(geo_accession = array_sample@header$geo_accession, 
                      sample_name = array_sample@header$source_name_ch1, 
                      sample_title = array_sample@header$title)
  
}) %>% 
  bind_rows(.)

# get RMA values
rma_tb <- purrr::map(names(array_list), function(array_name){
  
  # get one sample
  array_sample <- array_list[[array_name]] 
  
  # get RMA values
  rma_tb <- 
    array_sample %>% 
    Table(.) %>% 
    as_tibble(.) %>% 
    dplyr::mutate(geo_accession = array_sample@header$geo_accession)
  
}) %>% 
  bind_rows(.) %>% 
  tidyr::pivot_wider(id_cols = ID_REF, values_from = VALUE, names_from = geo_accession)


### PCA plot
# matrix for PCA
rma_matrix <- 
  rma_tb %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "ID_REF") %>% 
  as.matrix(.)

# calculates pca
pca <-
  rma_matrix %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes table for ggplot
pca_tb <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         geo_accession = colnames(rma_matrix)) %>%
  dplyr::left_join(sample_tb , by = "geo_accession")

### plot
# create bare plot object
pca_plot <- 
  ggplot() + 
  geom_point(data = pca_tb, aes(x = PC1, y = PC2, color = sample_name, fill = sample_name), size = 5, shape = 21) +
  # geom_label_repel(data = pca_tb, aes(x = PC1, y = PC2, label = sample_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5))) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "bottom")

# save labeled plot
ggsave(filename = file.path(outpath, "Barragan_2020_unpub_GSE152525.RMA.PCA_plot.png"),
       plot = pca_plot, width = 12, height = 10)


