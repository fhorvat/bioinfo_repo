### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Wed Nov 28 21:54:06 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis/smallRNA_clusters")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get PCA and variance of each PC
df2PCA <- function(df, sample_table, plot_name){
  
  # calculate PCA
  pca <-
    df %>%
    t(.) %>%
    stats::prcomp(.)
  
  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
  
  # makes data.frame for ggplot
  plot_df <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           sample_id = colnames(df)) %>%
    dplyr::left_join(sample_table , by = "sample_id") %>%
    dplyr::mutate(sample_id = str_remove_all(sample_id, "^s_|r") %>% str_replace_all(., "_", " "))
  
  # plot 
  pca_plot <- 
    ggplot(data = plot_df, aes(x = PC1, y = PC2, label = sample_id, fill = genotype, shape = transfection)) +
    geom_point(aes(fill = genotype), color = "black", size = 7.5) +
    scale_shape_manual(values  = c(21, 22)) +
    xlab(str_c("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(str_c("PC2: ", round(percentVar[2] * 100), "% variance")) +
    # guides(fill = FALSE, shape = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 23, size = 5)),
           shape = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = - 0.2),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, str_c("PCAplot.", plot_name, ".no_labels.png")), width = 12, height = 10)
  
  # add labels and save
  pca_plot <- 
    pca_plot + 
    ggrepel::geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", 
                              box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    ggsave(filename = file.path(outpath, str_c("PCAplot.", plot_name, ".with_labels.png")), width = 12, height = 10)
  
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# cluster expression path
cluster_path <- 
  list.files(inpath, pattern = "clusters.*\\.csv", full.names = T) %>% 
  .[!str_detect(., "mean")]

# sample table path
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Documentation/Eliska_mESC_MosIR.sampleTable.csv"

######################################################## READ DATA
# read cluster expression
cluster_df_list <- purrr::map(cluster_path, readr::read_csv)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
# mutate sample table
sample_table %<>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE$"))

### plot PCA based on expression of clusters
# loop through different clusters
for(n in 1:length(cluster_path)){
  
  # comparison name
  sample_name <- 
    cluster_path[n] %>% 
    basename(.) %>% 
    stringr::str_extract("RS10|RSP")
  
  # cluster PCA
  cluster_df <- 
    cluster_df_list[[n]] %>% 
    dplyr::select(-c(coordinates:class)) %>% 
    dplyr::mutate_all(funs(log(. + 1))) %>% 
    df2PCA(df = ., sample_table = sample_table, plot_name = str_c("clusters.", sample_name)) 
  
}

