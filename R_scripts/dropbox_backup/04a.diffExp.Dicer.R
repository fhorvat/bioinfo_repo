### INFO: 
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Stein_2015_PLoSGenet_GSE57514/Analysis")

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

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(plotly)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### experiment
# set experiment name
experiment <- "Stein_2015_PLoSGenet_GSE57514"

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10_new")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis")


### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.se.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.avgFPKM.csv$"), full.names = T)


### genome
# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# summarizedExperiment 
se <- readRDS(se_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)
  
######################################################## MAIN CODE
### prepare data
# filter ensembl genes info
genes_info_tidy <- 
  genes_info %>% 
  dplyr::mutate(cooridnates = str_c(seqnames, ":", start, "-", end, "|", strand)) %>% 
  dplyr::select(-c(seqnames:strand))

# get gene_id of protein coding genes
protein_genes <- 
  genes_info_tidy %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id


### exploratory analysis 
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  # dplyr::filter(str_detect(genotype, "Dicer")) %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
colnames(se_filt) <- str_remove(colnames(se_filt), ".genome.Aligned.sortedByCoord.out.bam|.total.bam")
se_filt <- se_filt[, colnames(se_filt) %in% sample_table_dds$sample_id]
colData(se_filt) <- DataFrame(sample_table_dds)

### DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~genotype)


### PCA plot
# data for PCA = rlog transformed counts
rlog_df <-
  rlog(dds, blind = T) %>%
  assay(.)

# calculates pca
pca <-
  rlog_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes data.frame for ggplot, plots PCA
pca_plot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " ")) %>%
  ggplot(data = ., aes(x = PC1, y = PC2, label = sample_id, fill = genotype)) +
  geom_point(aes(fill = genotype), color = "black", size = 7.5, shape = 21) +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  # scale_shape_manual(values  = c(21, 22, 24)) +
  # scale_fill_manual(values = c("red3", "blue3")) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  # guides(fill = FALSE, shape = FALSE) +
  # guides(shape = guide_legend(override.aes = list(size = 5))) +
  guides(fill = guide_legend(override.aes = list(shape = 23, size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, "results", str_c("PCAplot.", experiment, ".PC1_2.rlog.png")),
       plot = pca_plot, width = 10, height = 10)


### 3D PCA plot
# data for plot
pca_3Dplot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         PC3 = pca$x[, 3],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))

# make interactive 3D PCA plot
p <-
  plotly::plot_ly(data = pca_3Dplot,
                  x = ~PC1,
                  y = ~PC2,
                  z = ~PC3,
                  color = ~genotype
                  # colors = c("#1a75ff", "red3")
                  ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = str_c("PC1: ", round(percentVar[1] * 100), "% variance"), range = c(-40, 40)),
                      yaxis = list(title = str_c("PC2: ", round(percentVar[2] * 100), "% variance"), range = c(-40, 40)), 
                      zaxis = list(title = str_c("PC3: ", round(percentVar[3] * 100), "% variance"), range = c(-40, 40))))

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "results",  str_c("PCAplot.3D", experiment, ".PC1_2_3.rlog.html")),
                        selfcontained = T)


### distance heatmap
# calculate distance
dist_df <-
  rlog_df %>%
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
                   filename = file.path(outpath, "results", str_c("distHeatmap.", experiment, ".rlog.png")),
                   height = 10,
                   width = 12)


### DESeq2
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  dplyr::filter(str_detect(genotype, "Dicer")) %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
colnames(se_filt) <- str_remove(colnames(se_filt), ".genome.Aligned.sortedByCoord.out.bam|.total.bam")
se_filt <- se_filt[, colnames(se_filt) %in% sample_table_dds$sample_id]
colData(se_filt) <- DataFrame(sample_table_dds)

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~genotype)

# run DESeq
dds_deseq <- DESeq(dds)

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, contrast = c("genotype", "Dicer_KO", "Dicer_WT"))

# get results table
results_df <-
  dds_shrink %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  dplyr::arrange(padj) %>%
  dplyr::left_join(fpkm_tb %>% dplyr::select(gene_id, contains("Dicer")), by = "gene_id") %>%
  dplyr::left_join(genes_info_tidy, by = "gene_id") %>%
  dplyr::mutate(comparison = "Dicer.KO_vs_WT") %T>% 
  write_csv(., path = file.path(outpath, "results", str_c("diffExp", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "all.csv", sep = ".")))

# write only significant results, padj < 0.1
results_df %>%
  dplyr::filter(padj < 0.1) %T>%
  write_csv(., path = file.path(outpath, "results", str_c("diffExp", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "signif.csv", sep = ".")))


### DESeq2 MA plot
# data for plot
plot_df <- 
  results_df %>% 
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                sign = ifelse(padj < 0.1, "yes", "no"),
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, sign == "no", "not_sign"),
                regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>%
  dplyr::arrange(regulation)

# plot
ma_plot <- 
  ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
  geom_point(size = 3, shape = 20) +
  scale_x_log10(limits = c(1e-01, 1e5), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-7, 5),
                     breaks = c(-7:5)) +
  scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  guides(color = FALSE) +
  xlab("average expression") +
  ylab("log2FC") +
  ggtitle(str_c("Dicer", " KO vs. WT")) +
  theme_bw() +
  # theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
  #       axis.title.y = element_text(size = 15, vjust = 0.3)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, "results", str_c("MAplot", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)


