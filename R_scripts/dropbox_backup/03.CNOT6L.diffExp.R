### INFO: CNOT6L differential expression analysis
### DATE: Sun Mar 11 04:54:44 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis/test")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)

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

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
ensembl_tx_info <- readr::read_csv(file.path(inpath, "ensembl.Mus_musculus.GRCm38.89.20180305.UCSCnames.clean.transcriptInfo.csv"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(inpath, "ensembl.GRCm38.89.tx.CNOT6L.summarizedOverlaps.RDS"))

# read FPKM table
fpkm_df <- readr::read_csv(file = file.path(inpath, "ensembl.GRCm38.89.tx.CNOT6L.avgFPKM.csv"))

# sample table path
sample_table <- readr::read_csv(file = file.path(inpath, "CNOT6L.sample_table.csv"))

######################################################## MAIN CODE
### prepare data
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  dplyr::mutate(group = str_c(stage, "_", genotype)) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# get transcript_id of protein coding transcripts
protein_genes <- 
  ensembl_tx_info %>% 
  dplyr::filter(transcript_biotype == "protein_coding") %$%
  transcript_id

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
colData(se_filt) <- DataFrame(sample_table_dds)
se_filt$genotype <- factor(se_filt$genotype, levels = c("WT", "KO"))
se_filt$group <- factor(se_filt$group, levels = c("GV_WT", "GV_KO", "MII_WT", "MII_KO", "1C_WT", "1C_KO"))

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~group)

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
  ggplot(data = ., aes(x = PC1, y = PC2, label = sample_id, fill = genotype, shape = stage)) +
  geom_point(aes(fill = genotype), color = "black", size = 7.5) +
  # geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  scale_shape_manual(values  = c(21, 22, 24)) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  guides(fill = guide_legend(override.aes = list(shape = 23, size = 5)),
         shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, "results", str_c("PCAplot.CNOT6L.GRCm38.89.tx.rlog.png")),
       plot = pca_plot, width = 10, height = 10)

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
                   filename = file.path(outpath, "results", str_c("distHeatmap.CNOT6L.GRCm38.89.tx.rlog.png")),
                   height = 10,
                   width = 12)


### DESeq2
# run DESeq
dds_deseq <- DESeq(dds)

### get results
# loop through stages
for(filt_stage in unique(sample_table_dds$stage)){
  
  # compose groups
  groups <- c(str_c(filt_stage, "_KO"), str_c(filt_stage, "_WT"))
  
  # get results, shrink logFC
  dds_shrink <- lfcShrink(dds_deseq, contrast = c("group", groups))
  
  # get results table
  results_df <-
    dds_shrink %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(., var = "transcript_id") %>%
    as.tibble(.) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(comparison = str_c(filt_stage, ".KO_vs_WT")) %T>%
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.tx.all.csv")))
  
  # write only significant results, padj < 0.1
  results_df %>%
    dplyr::filter(padj < 0.1) %>%
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.tx.signif.csv")))
  
  ### MA plot
  # # read results 
  # results_df <- read_csv(file = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv")))
  
  # data for plot
  plot_df <- 
    results_df %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  sign = ifelse(padj < 0.1, "yes", "no"),
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, sign == "no", "not_sign"))
  
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
    # xlab("average expression") +
    # ylab("log2FC") +
    ggtitle(str_c(filt_stage, " KO vs. WT")) +
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
  ggsave(filename = file.path(outpath, "results", str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.tx.png")),
         plot = ma_plot, width = 10, height = 10)
  
}
