### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Green_2018_DevCell_GSE112393/Analysis/expression")

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

library(Seurat)
library(patchwork)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- file.path(getwd(), "Seurat_analysis")

# raw counts table
count_tb_path <- file.path(inpath, "GSE112393_MergedAdultMouseST25_DGE.filtered.matrix.RDS")

######################################################## READ DATA
# read filtered count data
count_tb <- readRDS(count_tb_path)

######################################################## MAIN CODE
### Seurat analysis
# create Seurat object
count_seurat <- CreateSeuratObject(counts = count_tb, project = "mouse_testis_scRNAseq")

# normalize data
count_norm <- NormalizeData(count_seurat, normalization.method = "LogNormalize", scale.factor = 10000)


### explore dataset
# find most variable genes
count_norm <- FindVariableFeatures(count_norm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(count_norm), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(count_norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(filename = file.path(outpath, "plot_01.most_variable_features.png"), width = 1500, height = 1000)
plot1 + plot2
dev.off()

# scale the data
all.genes <- rownames(count_norm)
count_norm <- ScaleData(count_norm, features = all.genes)

# perform linear dimensional reduction
count_norm <- RunPCA(count_norm, features = VariableFeatures(object = count_norm))

# plot PCA
png(filename = file.path(outpath, "plot_02.PCA.png"), width = 1000, height = 1000)
DimPlot(count_norm, reduction = "pca")
dev.off()

# plot heatmap
png(filename = file.path(outpath, "plot_03.heatmap.png"), width = 1000, height = 1000)
DimHeatmap(count_norm, dims = 1, cells = 500, balanced = TRUE)
dev.off()

# plot elbow plot
png(filename = file.path(outpath, "plot_04.elbow_PCA.png"), width = 1000, height = 1000)
ElbowPlot(count_norm)
dev.off()


### Clustering 
# run UMAP and tSNE
count_norm <- FindNeighbors(count_norm, dims = 1:10)
count_norm <- FindClusters(count_norm, resolution = 0.5)
count_norm <- RunUMAP(count_norm, dims = 1:10)
count_norm <- RunTSNE(count_norm, dims = 1:10)

# plot UMAP
png(filename = file.path(outpath, "plot_05.UMAP.png"), width = 1000, height = 1000)
DimPlot(count_norm, reduction = "umap")
dev.off()

# plot tSNE
png(filename = file.path(outpath, "plot_06.tSNE.png"), width = 1000, height = 1000)
DimPlot(count_norm, reduction = "tsne")
dev.off()


### Explore clusters
# find all markers of cluster 1
cluster1.markers <- FindMarkers(count_norm, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(count_norm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# plot violin plot - expression probability distributions across clusters
png(filename = file.path(outpath, "plot_07.violin.png"), width = 1000, height = 1000)
VlnPlot(count_norm, features = c("Insl3", "Cyp17a1", "Klk1", "Lcn2", "Cyp11a1", "Ptgds", "Lyz2", "Hsd3b1", "Klk1b21", "Hsd3b6"))
dev.off()

# plot violin plot - expression probability distributions across clusters
png(filename = file.path(outpath, "plot_08.tSNE_expression.png"), width = 1000, height = 1000)
FeaturePlot(count_norm, features = c("Insl3", "Cyp17a1", "Klk1", "Lcn2", "Cyp11a1", "Ptgds", "Lyz2", "Hsd3b1", "Klk1b21"), reduction = "tsne")
dev.off()

# expression heatmap
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
png(filename = file.path(outpath, "plot_09.expression_heatmap.png"), width = 1000, height = 1000)
DoHeatmap(count_norm, features = top10$gene) + NoLegend()
dev.off()