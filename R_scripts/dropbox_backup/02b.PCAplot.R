### INFO: Expression analysis mESC and oocytes sequenced in February 2018 and June 2018
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Analysis/expression")

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
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "mutate_cond.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Documentation/mESC_oocytes_2018.sample_table.csv"

# genes info
genes_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.geneInfo.csv"

######################################################## READ DATA
# read sample table
sample_table <- 
  readr::read_csv(file = sample_path) %>% 
  dplyr::filter(!(sample_id %in% c("s_ESC_DX_i3_JM7D1.SE", "s_MII_B6_MT_4.SE", "s_ESC_DXII_i6.SE"))) %>% 
  mutate_cond(stage == "ESC", genotype = str_remove_all(genotype, "\\.?[1,2]{1}$")) %>% 
  dplyr::mutate(genotype = str_replace(genotype, "i.*", "i"))

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read FPKM 
fpkm_df <- 
  readr::read_csv(file.path(outpath, "GRCm38.91.reducedExons.FPKM.csv")) %>% 
  dplyr::select(-c(seqnames:gene_description)) %>% 
  dplyr::select_at(.vars = vars(-matches("s_ESC_DX_i3_JM7D1.SE|s_MII_B6_MT_4.SE|s_ESC_DXII_i6.SE")))

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter FPKM to only protein coding genes
fpkm_df %<>% 
  dplyr::filter(gene_id %in% protein_genes)

# read summarizedExperiment from RDS file, filter to include only protein coding genes
se <- 
  readRDS(file = file.path(outpath, "GRCm38.91.reducedExons.summarizedOverlaps.RDS")) %>% 
  .[, -which(str_detect(colnames(.), "s_ESC_DX_i3_JM7D1.SE|s_MII_B6_MT_4.SE|s_ESC_DXII_i6.SE"))] %>% 
  .[rownames(.) %in% protein_genes, ]


### PCA plot - stages rlog
for(filt_stage in c("ESC|GV|MII", "GV|MII")){
  
  # filter sample table
  sample_table_filt <- 
    sample_table %>% 
    dplyr::filter(str_detect(stage, pattern = filt_stage)) %>% 
    as.data.frame(.) %>% 
    set_rownames(., .$short_name) 
  
  # filter summarizedExperiment to include only GV, set colData
  se_filt <- se[, colnames(se)[match(basename(sample_table_filt$bam_path), colnames(se))]]
  colData(se_filt) <- DataFrame(sample_table_filt)
  
  # make DESeqDataSet
  dds <- DESeqDataSet(se_filt, design = ~stage)
  
  # data for PCA
  PCA_df <-
    rlog(dds, blind = T) %>%
    # vst(dds, blind = T) %>%
    assay(.)
  
  # calculates pca
  pca <-
    PCA_df %>%
    t(.) %>%
    stats::prcomp(.)
  
  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
  
  # makes data.frame for ggplot, plots PCA
  PCA_plot <- 
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           short_name = colnames(PCA_df)) %>%
    dplyr::left_join(sample_table_filt, by = "short_name") %>%
    ggplot(data = ., aes(x = PC1, y = PC2, label = short_name, color = stage, shape = stage)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = short_name),
                    fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    scale_color_discrete(drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size = 5))) +
    xlab(str_c("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(str_c("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2)) +
    theme(axis.title.y = element_text(size = 15, vjust = 0.3))
  
  ggsave(filename = file.path(outpath, str_c("PCA.", str_replace_all(filt_stage, "\\|", "_"), ".rlog.stage.png")), 
         plot = PCA_plot, width = 10, height = 10)
  
  ### distance heatmap
  # calculate distance
  dist_df <-
    PCA_df %>%
    t(.) %>%
    dist(.)
  
  # make matrix
  dist_matrix <- as.matrix(dist_df)
  
  # set annotation column 
  anno_col <- sample_table_filt[, c("stage", "genotype"), drop = F]
  anno_colors <- list(stage = c(ESC = "#F8766D", GV = "#00BA38", MII = "#619CFF"))
  
  # plot
  pheatmap::pheatmap(dist_matrix,
                     cluster_rows = T, 
                     cluster_cols = T, 
                     clustering_distance_rows = dist_df,
                     clustering_distance_cols = dist_df,
                     # cellwidth = 15, 
                     # cellheight = 15, 
                     # col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                     # annotation_row = anno_col, 
                     annotation_colors = anno_colors,
                     filename = file.path(outpath, str_c("dist.", str_replace_all(filt_stage, "\\|", "_"), ".rlog.stage.png")),
                     height = 10,
                     width = 10)
  
}


### PCA plot - genotype rlog
for(filt_stage in c("ESC", "GV", "MII")){
  
  filt_stage <- "ESC"
  
  # filter sample table
  sample_table_filt <- 
    sample_table %>% 
    dplyr::filter(str_detect(stage, pattern = filt_stage))  %>% 
    as.data.frame(.) %>% 
    set_rownames(., .$short_name)  
  
  # filter summarizedExperiment to include only GV, set colData
  se_filt <- se[, colnames(se)[match(basename(sample_table_filt$bam_path), colnames(se))]]
  colData(se_filt) <- DataFrame(sample_table_filt)
  
  # make DESeqDataSet
  dds <- DESeqDataSet(se_filt, design = ~genotype)
  
  # data for PCA
  PCA_df <-
    rlog(dds, blind = T) %>%
    assay(.)
  
  # calculates pca
  pca <-
    PCA_df %>%
    t(.) %>%
    stats::prcomp(.)
  
  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
  
  # makes data.frame for ggplot, plots PCA
  PCA_plot <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           short_name = colnames(PCA_df)) %>%
    dplyr::left_join(sample_table_filt, by = "short_name") %>%
    ggplot(data = ., aes(x = PC1, y = PC2, label = short_name, color = genotype, shape = genotype)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = short_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    scale_colour_manual(values = c("black", "red3", "#1a75ff", "#ee34ff", "orange", "purple", "cyan", "darkgreen")) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 2, 5)) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2)) +
    theme(axis.title.y = element_text(size = 15, vjust = 0.3))
  
  ggsave(filename = file.path(outpath, str_c("PCA.", filt_stage, ".rlog.genotype.png")),
         plot = PCA_plot, width = 10, height = 10)
  
  ### distance heatmap
  # calculate distance
  dist_df <-
    PCA_df %>%
    t(.) %>%
    dist(.)
  
  # make matrix
  dist_matrix <- as.matrix(dist_df)
  colnames(dist_matrix) <- NULL
  
  # set annotation column 
  anno_col <- sample_table_filt[, "genotype", drop = F]
  
  # plot
  pheatmap::pheatmap(dist_matrix,
                     cluster_rows = T, 
                     cluster_cols = T, 
                     clustering_distance_rows = dist_df,
                     clustering_distance_cols = dist_df,
                     # cellwidth = 15, 
                     # cellheight = 15, 
                     col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                     # annotation_row = anno_col, 
                     # annotation_colors = anno_colors,
                     filename = file.path(outpath, str_c("dist.", filt_stage, ".rlog.genotype.png")),
                     height = 10,
                     width = 10)
  
}


### PCA plot - stages log2FPKM
for(filt_stage in c("ESC|GV|MII", "GV|MII")){

  filt_stage <- "ESC|GV|MII"
  
  # filter sample table
  sample_table_filt <-
    sample_table %>%
    dplyr::filter(str_detect(stage, pattern = filt_stage))

  # set vector for renaming columns of FPKM data.frame
  column_names <-
    sample_table_filt$short_name %>%
    magrittr::set_names(., sample_table_filt$sample_id)

  # filter summarizedExperiment to include only GV, set colData
  PCA_df <-
    fpkm_df %>%
    dplyr::select_at(.vars = vars(matches(filt_stage))) %>%
    dplyr::mutate_all(.funs = funs(log2(. + 1))) %>%
    magrittr::set_colnames(column_names[colnames(.)])

  # calculates pca
  pca <-
    PCA_df %>%
    t(.) %>%
    stats::prcomp(.)

  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

  # makes data.frame for ggplot, plots PCA
  PCA_plot <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           short_name = colnames(PCA_df)) %>%
    dplyr::left_join(sample_table_filt, by = "short_name") %>%
    ggplot(data = ., aes(x = PC1, y = PC2, label = short_name, color = stage, shape = stage)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = short_name),
                    fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    scale_color_discrete(drop = FALSE) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2)) +
    theme(axis.title.y = element_text(size = 15, vjust = 0.3))

  ggsave(filename = file.path(outpath, str_c("PCA.", str_replace_all(filt_stage, "\\|", "_"), ".logFPKM.stage.png")),
         plot = PCA_plot, width = 10, height = 10)

  ### distance heatmap
  # calculate distance
  dist_df <-
    PCA_df %>%
    t(.) %>%
    dist(.)

  # make matrix
  dist_matrix <- as.matrix(dist_df)

  # set annotation column
  anno_col <- sample_table_filt[, c("stage", "genotype"), drop = F]
  anno_colors <- list(stage = c(ESC = "#F8766D", GV = "#00BA38", MII = "#619CFF"))

  # plot
  pheatmap::pheatmap(dist_matrix,
                     cluster_rows = T,
                     cluster_cols = T,
                     clustering_distance_rows = dist_df,
                     clustering_distance_cols = dist_df,
                     # cellwidth = 15,
                     # cellheight = 15,
                     col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                     # annotation_row = anno_col,
                     # annotation_colors = anno_colors,
                     filename = file.path(outpath, str_c("dist.", str_replace_all(filt_stage, "\\|", "_"), ".logFPKM.stage.png")),
                     height = 10,
                     width = 10)

}


### PCA plot - genotype log2FPKM
for(filt_stage in c("GV", "MII", "ESC")){
  
  filt_stage <- "ESC"
  
  # filter sample table
  sample_table_filt <-
    sample_table %>%
    dplyr::filter(str_detect(stage, pattern = filt_stage))

  # set vector for renaming columns of FPKM data.frame
  column_names <-
    sample_table_filt$short_name %>%
    magrittr::set_names(., sample_table_filt$sample_id)

  # filter summarizedExperiment to include only GV, set colData
  PCA_df <-
    fpkm_df %>%
    dplyr::select_at(.vars = vars(matches(filt_stage))) %>%
    dplyr::mutate_all(.funs = funs(log2(. + 1))) %>%
    magrittr::set_colnames(column_names[colnames(.)])

  # calculates pca
  pca <-
    PCA_df %>%
    t(.) %>%
    stats::prcomp(.)

  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

  # makes data.frame for ggplot, plots PCA
  PCA_plot <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           short_name = colnames(PCA_df)) %>%
    dplyr::left_join(sample_table_filt, by = "short_name") %>%
    ggplot(data = ., aes(x = PC1, y = PC2, label = short_name, color = genotype, shape = genotype)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = short_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    scale_colour_manual(values = c("black", "red3", "#1a75ff", "#ee34ff", "orange")) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2)) +
    theme(axis.title.y = element_text(size = 15, vjust = 0.3))

  ggsave(filename = file.path(outpath, str_c("PCA.", str_replace_all(filt_stage, "\\|", "_"), ".logFPKM.genotype.png")),
         plot = PCA_plot, width = 10, height = 10)

  ### distance heatmap
  # calculate distance
  dist_df <-
    PCA_df %>%
    t(.) %>%
    dist(.)

  # make matrix
  dist_matrix <- as.matrix(dist_df)

  # set annotation column
  anno_col <- sample_table_filt[, c("stage", "genotype"), drop = F]
  anno_colors <- list(stage = c(ESC = "#F8766D", GV = "#00BA38", MII = "#619CFF"))

  # plot
  pheatmap::pheatmap(dist_matrix,
                     cluster_rows = T,
                     cluster_cols = T,
                     clustering_distance_rows = dist_df,
                     clustering_distance_cols = dist_df,
                     # cellwidth = 15,
                     # cellheight = 15,
                     col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                     # annotation_row = anno_col,
                     # annotation_colors = anno_colors,
                     filename = file.path(outpath, str_c("dist.", str_replace_all(filt_stage, "\\|", "_"), ".logFPKM.genotype.png")),
                     height = 10,
                     width = 10)

}


### PCA plot - new mESC samples sequenced in June 2018
# filter sample table
sample_table_filt <- 
  sample_table %>% 
  dplyr::filter(stringr::str_detect(sample_id, "DXII"))  %>% 
  as.data.frame(.) %>% 
  set_rownames(., .$short_name)  

# filter summarizedExperiment to include only GV, set colData
se_filt <- se[, colnames(se)[match(basename(sample_table_filt$bam_path), colnames(se))]]
colData(se_filt) <- DataFrame(sample_table_filt)

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~genotype)

# data for PCA
PCA_df <-
  rlog(dds, blind = T) %>%
  assay(.)

# calculates pca
pca <-
  PCA_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes data.frame for ggplot, plots PCA
PCA_plot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         short_name = colnames(PCA_df)) %>%
  dplyr::left_join(sample_table_filt, by = "short_name") %>%
  ggplot(data = ., aes(x = PC1, y = PC2, label = short_name, color = genotype, shape = genotype)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = short_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  scale_colour_manual(values = c("black", "red3", "#1a75ff", "#ee34ff")) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2)) +
  theme(axis.title.y = element_text(size = 15, vjust = 0.3))

ggsave(filename = file.path(outpath, "PCA.mESC.June2018.rlog.genotype.png"),
       plot = PCA_plot, width = 10, height = 10)

### distance heatmap
# calculate distance
dist_df <-
  PCA_df %>%
  t(.) %>%
  dist(.)

# make matrix
dist_matrix <- as.matrix(dist_df)
colnames(dist_matrix) <- NULL

# set annotation column 
anno_col <- sample_table_filt[, "genotype", drop = F]

# plot
pheatmap::pheatmap(dist_matrix,
                   cluster_rows = T, 
                   cluster_cols = T, 
                   clustering_distance_rows = dist_df,
                   clustering_distance_cols = dist_df,
                   # cellwidth = 15, 
                   # cellheight = 15, 
                   col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   # annotation_row = anno_col, 
                   # annotation_colors = anno_colors,
                   filename = file.path(outpath, "dist.mESC.June2018.rlog.genotype.png"),
                   height = 10,
                   width = 10)