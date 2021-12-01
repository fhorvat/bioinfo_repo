### INFO: PCA plot
### DATE: Tue Dec 11 22:50:27 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(ggrepel)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

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
    ggplot(data = plot_df, aes(x = PC1, y = PC2, label = sample_id, fill = genotype, shape = stage)) +
    geom_point(aes(fill = genotype), color = "black", size = 7.5) +
    scale_shape_manual(values  = c(21, 22, 24)) +
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
    ggsave(filename = file.path(outpath, str_c("PCAplot.no_labels.", plot_name, ".png")), width = 12, height = 10)
  
  # add labels and save
  pca_plot <- 
    pca_plot + 
    ggrepel::geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", 
                              box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    ggsave(filename = file.path(outpath, str_c("PCAplot.with_labels.", plot_name, ".png")), width = 12, height = 10)
  
  
}

######################################################## PATH VARIABLES
# set input path
inpath <- getwd()

# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation/lnc1_KO.RNAseq.20181211.sampleTable.clean.csv"

# genes info
genes_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.geneInfo.csv"

# summarizedExperiment path
se_path <- file.path(inpath, "Lnc1_KO.2018_Dec.GRCm38.91.reducedExons.summarizedOverlaps.RDS")

# FPKM path
fpkm_path <- file.path(inpath, "Lnc1_KO.2018_Dec.GRCm38.91.reducedExons.FPKM.csv")

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path) 

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read summarizedExperiment
se <- readRDS(file = se_path) 

# read FPKM 
fpkm_df <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
# filter sample table
sample_table_dds <- 
  sample_table %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE")) %>% 
  tidyr::unite(col = "group", genotype, stage, sep = "_", remove = F) %>% 
  as.data.frame(.) %>% 
  set_rownames(., .$sample_id) 

# # get gene_id of protein coding genes
# protein_genes <- 
#   genes_info %>% 
#   dplyr::filter(gene_biotype == "protein_coding") %$%
#   gene_id


### summarizedOverlaps 
# change colnames
colnames(se) <- str_remove(colnames(se), "\\.SE\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam")

# read summarizedExperiment from RDS file, filter to include only protein coding genes
se_filt <- 
  se %>% 
  # .[rownames(.) %in% protein_genes, ] %>% 
  .[, colnames(.)[match(rownames(sample_table_dds), colnames(.))]]

# add column data to SE
colData(se_filt) <- DataFrame(sample_table_dds)

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~group) 

# data for PCA = rlog transformed counts
pca_plot <-
  rlog(dds, blind = T) %>%
  assay(.)

# calculate PCA and plot
df2PCA(df = pca_plot, sample_table = sample_table_dds, plot_name = "Lnc1_KO.2018_Dec.all_biotypes.rlog") 

# calculate distance
dist_df <-
  pca_plot %>%
  t(.) %>%
  dist(.)

# make matrix
dist_matrix <- as.matrix(dist_df)
colnames(dist_matrix) <- NULL

# plot
pheatmap::pheatmap(dist_matrix,
                   clustering_distance_rows = dist_df,
                   clustering_distance_cols = dist_df,
                   col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   filename = file.path(outpath, "distHeatmap.Lnc1_KO.2018_Dec.all_biotypes.rlog.png"),
                   height = 10,
                   width = 12)

### FPKM
# filter FPKM
pca_plot <- 
  fpkm_df %>% 
  # dplyr::filter(gene_id %in% protein_genes) %>% 
  dplyr::select(-c(gene_id, coordinates:gene_description)) %>% 
  dplyr::mutate_all(funs(log(. + 1)))

# calculate PCA and plot
df2PCA(df = pca_plot, sample_table = sample_table_dds, plot_name = "Lnc1_KO.2018_Dec.all_biotypes.FPKM")

# calculate distance
dist_df <-
  pca_plot %>%
  t(.) %>%
  dist(.)

# make matrix
dist_matrix <- as.matrix(dist_df)
colnames(dist_matrix) <- NULL

# plot
pheatmap::pheatmap(dist_matrix,
                   clustering_distance_rows = dist_df,
                   clustering_distance_cols = dist_df,
                   col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   filename = file.path(outpath, "distHeatmap.Lnc1_KO.2018_Dec.all_biotypes.FPKM.png"),
                   height = 10,
                   width = 12)
