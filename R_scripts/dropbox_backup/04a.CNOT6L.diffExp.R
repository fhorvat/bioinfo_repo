### INFO: CNOT6L differential expression analysis
### DATE: Sun Mar 11 04:54:44 2018
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
library(plotly)

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

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

# # reduced exons path
# exons_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.reducedExons.RDS"

######################################################## READ DATA
# # read ENSEMBL reduced exons
# exons_gr <- readRDS(file = exons_path)

# read info about all ENSEMBL annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(inpath, "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"))

# read FPKM table
fpkm_df <- readr::read_csv(file = file.path(inpath, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

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

# filter ensembl genes info
ensembl_genes_info_filt <- 
  ensembl_genes_info %>% 
  dplyr::mutate(cooridnates = str_c(seqnames, ":", start, "-", end, "|", strand)) %>% 
  dplyr::select(-c(seqnames:strand))

# get gene_id of protein coding genes
protein_genes <- 
  ensembl_genes_info_filt %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# # get total length of all exons for each transcript
# exons_width <-
#   width(exons_gr) %>%
#   sum(.) %>%
#   tibble(gene_id = names(.), width = .)

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
  # scale_fill_manual(values = c("red3", "blue3")) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  # guides(fill = FALSE, shape = FALSE) +
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
ggsave(filename = file.path(outpath, "results", str_c("PCAplot.CNOT6L.GRCm38.89.PC1_2.rlog.png")),
       plot = pca_plot, width = 10, height = 10)

### 3D PCA plot
# data for plot
pca_3Dplot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         PC3 = pca$x[, 3],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "), 
                stage = factor(stage, levels = c("GV", "MII", "1C")), 
                genotype = factor(genotype, levels = c("WT", "KO")))

# make interactive 3D PCA plot
p <-
  plotly::plot_ly(data = pca_3Dplot,
                  x = ~PC1,
                  y = ~PC2,
                  z = ~PC3,
                  symbol = ~stage,
                  symbols = c("square", "diamond", "circle"),
                  color = ~genotype,
                  colors = c("#1a75ff", "red3")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = str_c("PC1: ", round(percentVar[1] * 100), "% variance"), range = c(-40, 40)),
                      yaxis = list(title = str_c("PC2: ", round(percentVar[2] * 100), "% variance"), range = c(-40, 40)), 
                      zaxis = list(title = str_c("PC3: ", round(percentVar[3] * 100), "% variance"), range = c(-40, 40))))

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(p),
                        file = file.path(outpath, "results", str_c("PCAplot.3D.CNOT6L.GRCm38.89.PC1_2_3.rlog.html")),
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
                   filename = file.path(outpath, "results", str_c("distHeatmap.CNOT6L.GRCm38.89.rlog.png")),
                   height = 10,
                   width = 12)


### DESeq2
# run DESeq
dds_deseq <- DESeq(dds)

### get results
# loop through stages
for(filt_stage in unique(sample_table_dds$stage)){
  
  # # compose groups
  # groups <- c(str_c(filt_stage, "_KO"), str_c(filt_stage, "_WT"))
  # 
  # # get results, shrink logFC
  # dds_shrink <- lfcShrink(dds_deseq, contrast = c("group", groups))
  # 
  # # get results table
  # results_df <-
  #   dds_shrink %>%
  #   as.data.frame(.) %>%
  #   tibble::rownames_to_column(., var = "gene_id") %>%
  #   as.tibble(.) %>%
  #   dplyr::arrange(padj) %>%
  #   dplyr::left_join(fpkm_df %>% dplyr::select(gene_id, matches(filt_stage)), by = "gene_id") %>%
  #   dplyr::left_join(ensembl_genes_info_filt, by = "gene_id") %>%
  #   dplyr::mutate(comparison = str_c(filt_stage, ".KO_vs_WT")) %T>%
  #   write_csv(., path = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv")))
  # 
  # # write only significant results, padj < 0.1
  # results_df %>%
  #   dplyr::filter(padj < 0.1) %T>%
  #   write_csv(., path = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.signif.csv")))
  
  
  ### DESeq2 MA plot
  # read results 
  results_df <- read_csv(file = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv")))
  
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
  ggsave(filename = file.path(outpath, "results", str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.png")),
         plot = ma_plot, width = 10, height = 10)
  
  
  # ### FPKM MA plot
  # # data.frame for plot
  # plot_df <- 
  #   results_df %>% 
  #   dplyr::select(matches(filt_stage), padj, gene_id) %>% 
  #   dplyr::rename_at(.vars = vars(matches("KO|WT")), .funs = funs(c("lfc_KO", "lfc_WT"))) %>% 
  #   dplyr::mutate(mean = dplyr::select(., lfc_KO, lfc_WT) %>% rowMeans(., na.rm = T), 
  #                 lfc = (log2(lfc_KO + 1) - log2(lfc_WT + 1)), 
  #                 padj = replace(padj, is.na(padj), 1), 
  #                 sign = ifelse(padj < 0.1, "yes", "no"),
  #                 regulation = ifelse(lfc > 0, "up", "down"), 
  #                 regulation = replace(regulation, sign == "no", "not_sign"), 
  #                 regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>%
  #   dplyr::arrange(regulation)
  # 
  # # plot
  # ma_plot <- 
  #   ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
  #   geom_point(size = 3, shape = 20) +
  #   scale_x_log10(limits = c(1e-03, 1e5), 
  #                 breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  #   scale_y_continuous(limits = c(-7, 5),
  #                      breaks = c(-7:5)) +
  #   scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  #   guides(color = FALSE) +
  #   # xlab("average expression") +
  #   # ylab("log2FC") +
  #   ggtitle(str_c(filt_stage, " KO vs. WT")) +
  #   theme_bw() +
  #   # theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
  #   #       axis.title.y = element_text(size = 15, vjust = 0.3)) +
  #   theme(axis.title.x = element_blank(), 
  #         axis.title.y = element_blank()) +
  #   theme(axis.text.x = element_text(size = 15), 
  #         axis.text.y = element_text(size = 15),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # # save plot
  # ggsave(filename = file.path(outpath, "results", str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT", ".FPKM.GRCm38.89.png")),
  #        plot = ma_plot, width = 10, height = 10)
  # 
  # 
  # ### FPM MA plot
  # fpm_df <-
  #   plot_df <- 
  #   results_df %>% 
  #   dplyr::select(matches(filt_stage), gene_id) %>% 
  #   dplyr::rename_at(.vars = vars(matches("KO|WT")), .funs = funs(c("lfc_KO", "lfc_WT"))) %>% 
  #   tidyr::gather(key = sample_id, value = fpkm, -gene_id) %>%
  #   dplyr::left_join(., exons_width, by = "gene_id") %>%
  #   dplyr::mutate(width = round(width / 1E3, 3), 
  #                 fpm = fpkm * width) %>%
  #   dplyr::select(gene_id, sample_id, fpm) %>%
  #   tidyr::spread(key = sample_id, value = fpm)
  # 
  # # data.frame for plot
  # plot_df <- 
  #   results_df %>% 
  #   dplyr::select(padj, gene_id) %>% 
  #   dplyr::left_join(., fpm_df, by = "gene_id") %>% 
  #   dplyr::mutate(mean = dplyr::select(., lfc_KO, lfc_WT) %>% rowMeans(., na.rm = T), 
  #                 lfc = (log2(lfc_KO + 1) - log2(lfc_WT + 1)), 
  #                 padj = replace(padj, is.na(padj), 1), 
  #                 sign = ifelse(padj < 0.1, "yes", "no"),
  #                 regulation = ifelse(lfc > 0, "up", "down"), 
  #                 regulation = replace(regulation, sign == "no", "not_sign"), 
  #                 regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>%
  #   dplyr::arrange(regulation)
  #   
  # # plot
  # ma_plot <- 
  #   ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) + 
  #   geom_point(size = 3, shape = 20) +
  #   scale_x_log10(limits = c(1e-03, 1e5), 
  #                 breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  #   scale_y_continuous(limits = c(-7, 5),
  #                      breaks = c(-7:5)) +
  #   scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  #   guides(color = FALSE) +
  #   # xlab("average expression") +
  #   # ylab("log2FC") +
  #   ggtitle(str_c(filt_stage, " KO vs. WT")) +
  #   theme_bw() +
  #   # theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
  #   #       axis.title.y = element_text(size = 15, vjust = 0.3)) +
  #   theme(axis.title.x = element_blank(), 
  #         axis.title.y = element_blank()) +
  #   theme(axis.text.x = element_text(size = 15), 
  #         axis.text.y = element_text(size = 15),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # # save plot
  # ggsave(filename = file.path(outpath, "results", str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT", ".FPM.GRCm38.89.png")),
  #        plot = ma_plot, width = 10, height = 10)
  
}



#### MII vs. GV WT
# compose groups
groups <- c("MII_WT", "GV_WT")

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, contrast = c("group", groups))

# get results table
plot_df <-
  dds_shrink %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  dplyr::arrange(padj) %>%
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
  ggtitle("MII vs GV WT") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, "results", "MAplot.CNOT6L.MII_WT.vs.GV_WT.GRCm38.89.png"),
       plot = ma_plot, width = 10, height = 10)