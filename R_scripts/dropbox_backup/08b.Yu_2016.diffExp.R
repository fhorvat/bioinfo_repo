### INFO: Yu 2016 BTG4 KO differential expression analysis
### DATE: Tue Apr 10 00:05:56 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis")

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

# set experiment name
experiment_name <- "Yu2016"

# ENSEMBL annotated genes info path
genes_info_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.geneInfo.csv"

# summarizedExperiment RDS path
se_path <- list.files(path = inpath, pattern = str_c(".*", experiment_name, ".se.RDS"))

# FPKM table path
fpkm_path <- list.files(path = inpath, pattern = str_c(".*", experiment_name, ".avgFPKM.csv"))

######################################################## READ DATA
# read info about all ENSEMBL annotated genes
genes_info <- readr::read_csv(genes_info_path)

# read summarizedExperiment from RDS file
se <- readRDS(file = se_path)

# read FPKM table
fpkm_df <- readr::read_csv(file = fpkm_path)

# sample table path
sample_table <- readr::read_csv(file = file.path(inpath, str_c(experiment_name, ".sampleTable.csv")))

######################################################## MAIN CODE
### prepare data
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# prepare sample table for dds
sample_table_dds <-
  sample_table %>%
  dplyr::mutate(sample_id = str_replace_all(sample_id, ".PE", ""),
                group = str_c(stage, "_", genotype)) %>%
  as.data.frame(.) %>%
  magrittr::set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
se_filt <- se_filt[, colnames(se_filt)[match(basename(sample_table_dds$bam_path), colnames(se_filt))]]
colData(se_filt) <- DataFrame(sample_table_dds)
se_filt$stage <- factor(se_filt$stage, levels = c("GV", "MII", "1C"))
se_filt$genotype <- factor(se_filt$genotype, levels = c("KO", "WT"))

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
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>%
                  str_replace_all(., "_", " ")) %>%
  ggplot(data = ., aes(x = PC1, y = PC2, label = sample_id, fill = genotype, shape = stage)) +
  geom_point(aes(fill = genotype), color = "black", size = 7.5) +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
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
ggsave(filename = file.path(outpath, "results", str_c("PCAplot.", experiment_name, ".GRCm38.89.rlog.png")),
       plot = pca_plot, width = 10, height = 10)


### DESeq2
# run DESeq
dds_deseq <- DESeq(dds)

### get results, plot MA and FPKM MA
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
    tibble::rownames_to_column(., var = "gene_id") %>%
    as.tibble(.) %>%
    dplyr::arrange(padj) %>%
    dplyr::left_join(fpkm_df %>% dplyr::select(gene_id, matches(filt_stage)), by = "gene_id") %>%
    dplyr::left_join(genes_info, by = "gene_id") %>%
    dplyr::mutate(comparison = str_c(filt_stage, ".KO_vs_WT")) %T>%
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.", experiment_name, ".", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv")))
  
  # # write only significant results, padj < 0.1
  # results_df %>%
  #   dplyr::filter(padj < 0.1) %T>%
  #   write_csv(., path = file.path(outpath, "results", str_c("diffExp.", experiment_name, ".", filt_stage, ".KO_vs_WT", ".GRCm38.89.signif.csv")))
  # 
  # 
  # ### MA plot
  # # data for plot
  # plot_df <- 
  #   results_df %>% 
  #   dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
  #   dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
  #                 sign = ifelse(padj < 0.1, "yes", "no"),
  #                 regulation = ifelse(lfc > 0, "up", "down"), 
  #                 regulation = replace(regulation, sign == "no", "not_sign"), 
  #                 regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>% 
  #   dplyr::arrange(regulation)
  # 
  # # ENSMUSG00000032056 (Btg4) data
  # plot_btg4 <- 
  #   plot_df %>% 
  #   dplyr::filter(gene_id == "ENSMUSG00000032056")
  # 
  # # plot
  # ma_plot <- 
  #   ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
  #   geom_point(size = 3, shape = 20) +
  #   geom_label(data = plot_btg4, aes(x = mean, y = lfc), label = "Btg4", hjust = 0, color = "black", nudge_x = 0.1) +
  #   geom_point(data = plot_btg4, aes(x = mean, y = lfc), size = 6, shape = 1, color = "black") + 
  #   scale_shape_discrete(solid = F) +
  #   scale_x_log10(limits = c(1e-01, 1e5), 
  #                 breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  #   scale_y_continuous(limits = c(-7, 5),
  #                      breaks = c(-7:5)) +
  #   scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  #   guides(color = FALSE) +
  #   ggtitle(str_c(filt_stage, " KO vs. WT")) +
  #   theme_bw() +
  #   theme(axis.title.x = element_blank(), 
  #         axis.title.y = element_blank()) +
  #   theme(axis.text.x = element_text(size = 15), 
  #         axis.text.y = element_text(size = 15),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # # save plot
  # ggsave(filename = file.path(outpath, "results", str_c("MAplot.", experiment_name, ".", filt_stage, ".KO_vs_WT", ".GRCm38.89.png")),
  #        plot = ma_plot, width = 10, height = 10)
  
  
  ### FPKM MA plot
  # data.frame for plot
  plot_df <- 
    results_df %>% 
    dplyr::select(matches(filt_stage), padj, gene_id, log2FoldChange) %>% 
    dplyr::rename_at(.vars = vars(matches("KO|WT")), .funs = funs(c("lfc_KO", "lfc_WT"))) %>% 
    dplyr::mutate(mean = dplyr::select(., lfc_KO, lfc_WT) %>% rowMeans(., na.rm = T), 
                  lfc = (log2(lfc_KO + 1) - log2(lfc_WT + 1)), 
                  padj = replace(padj, is.na(padj), 1), 
                  sign = ifelse(padj < 0.1, "yes", "no"),
                  regulation = ifelse(log2FoldChange > 0, "up", "down"), 
                  regulation = replace(regulation, sign == "no", "not_sign"), 
                  regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>%
    dplyr::arrange(regulation)
  
  # plot
  ma_plot <- 
    ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
    geom_point(size = 5, shape = 20) +
    geom_hline(yintercept = 0) + 
    scale_x_log10(limits = c(0.001, 1e5), 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_continuous(breaks = c(-5:3)) +
    scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
    guides(color = FALSE) +
    xlab("average expression") +
    ylab(str_c("log2(", filt_stage, " KO vs. WT)")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
          axis.title.y = element_text(size = 15, vjust = 0.3), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("MAplot.", experiment_name, ".FPKM.", filt_stage, ".KO_vs_WT", ".GRCm38.89.png")),
         plot = ma_plot, width = 10, height = 10)
  
}
