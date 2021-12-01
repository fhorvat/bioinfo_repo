### INFO:
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_mESC.Dicer_mutants.small_RNAseq.2021_May/Analysis/expression")

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
library(openxlsx)

library(viridis)
library(dendsort)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# sort clusters
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_mESC.Dicer_mutants.small_RNAseq.2021_May"

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# counts path
counts_path <- file.path(inpath, "miRBase.22.mm10.20181605.counts.txt")

# FPM path
fpm_mean_path <- file.path(inpath, "miRBase.22.mm10.20181605.FPM_mean.csv")
fpm_path <- file.path(inpath, "miRBase.22.mm10.20181605.FPM.csv")

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.))) %>%
  dplyr::mutate(Geneid = make.unique(Geneid))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read FPM tables
fpm_mean_tb <- readr::read_csv(fpm_mean_path)
fpm_tb <- readr::read_csv(fpm_path)

######################################################## MAIN CODE
## set results
results_groups <- c("")
results_groups <- str_split(results_groups, pattern = " ") %>% unlist()

# prepare sample table for DESeq colData
sample_table_dds <-
  sample_table %>%
  dplyr::filter(mutation == "null_mESCs") %>%
  dplyr::filter(sample_id %in% str_remove_all(colnames(counts_tb), "\\.24to31nt|\\.19to32nt|\\.21to23nt|\\.bam$")) %>% 
  dplyr::select(sample_id, mutation) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# counts table
se <-
  counts_tb %>%
  dplyr::select(-c(Chr:Length)) %>%
  dplyr::rename(gene_id = Geneid) %>%
  dplyr::mutate_if(is.numeric, round, digits = 0) %>%
  as.data.frame(.) %>%
  set_rownames(., .$gene_id) %>%
  dplyr::select(-gene_id) %>%
  as.matrix(.)

# filter summarizedExperiment to include only chosen stage, set colData
se_filt <- se
se_filt <- SummarizedExperiment(list(counts = se_filt))
colnames(se_filt) <- str_remove(colnames(se_filt), "\\.21to23nt\\.bam$")
se_filt <- se_filt[, colnames(se_filt)[match(rownames(sample_table_dds), colnames(se_filt))]]

# check if colnames of assay match rownames of sample table DataFrame
if(all(colnames(se_filt) == rownames(sample_table_dds))){
  
  # set sample table as colData
  colData(se_filt) <- DataFrame(sample_table_dds)
  
}else{
  
  # stop script with warrning
  stop("Columns in assay are not matching row of sample table. Please check your data annotation")
  
}

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~1)

# data for PCA = rlog transformed counts
rlog_df <-
  rlog(dds, blind = T) %>%
  assay(.)


### PCA plot
if(TRUE){
  
  # calculates pca
  pca <-
    rlog_df %>%
    t(.) %>%
    .[ , which(apply(., 2, var) != 0)] %>%
    stats::prcomp(., scale. = T)
  
  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
  
  # makes table for ggplot
  pca_tb <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           sample_id = colnames(rlog_df)) %>%
    dplyr::left_join(sample_table_dds , by = "sample_id") %>%
    dplyr::mutate(sample_id = str_remove_all(sample_id, "s_|r[0-9]+\\.SE|Dicer_SOM") %>% 
                    str_replace_all(., "_", " ") %>% 
                    str_trim(.)) %>% 
    dplyr::mutate(mutation = factor(mutation, levels = c("delA", "P197A", "P216T", 
                                                         "A237V", "E176Q", "K70N", 
                                                         "G69V", "delHelicase")))
  
  # set axis limits
  axis_lim <-
    c(pca_tb$PC1, pca_tb$PC2) %>%
    abs(.) %>%
    max(.) %>%
    ceiling(.)
  
  ### plot
  # create bare plot object
  pca_plot <- 
    ggplot(data = pca_tb, aes(x = PC1, y = PC2)) + 
    geom_point(aes(fill = mutation), color = "black", size = 6, shape = 21) +
    guides(fill = guide_legend(override.aes = list(shape = 23, size = 5))) + 
    scale_x_continuous(limits = c(-axis_lim, axis_lim)) +
    scale_y_continuous(limits = c(-axis_lim, axis_lim)) +
    scale_fill_manual(values = c("delA" = "#00B050", "P197A" = "#0000FF", "P216T" = "#FF0000", 
                                 "A237V" = "#FFC000", "E176Q" = "#00B0F0", "K70N" = "#A6A6A6", 
                                 "G69V" = "#C6E0B4", "delHelicase" = "#FF00FF")) +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13, angle = 90),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.title = element_blank()) +
    theme(legend.position = "right")
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("miRBase", "plot", "PCA.PC1_PC2", "rlog", "mutation", "png", sep = ".")),
         plot = pca_plot, width = 10, height = 10)
  
  # add labels
  pca_plot <- 
    pca_plot + 
    geom_label_repel(aes(label = sample_id),
                     fontface = "bold", color = "black", box.padding = 0.35,
                     point.padding = 0.5, segment.color = "grey50")
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("miRBase", "plot", "PCA.PC1_PC2", "rlog", "mutation", "labeled", "png", sep = ".")),
         plot = pca_plot, width = 10, height = 10)
  
}

### distance heatmap
if(TRUE){
  
  # calculate distance
  dist_df <-
    rlog_df %>%
    t(.) %>%
    dist(.)
  
  # make matrix
  dist_matrix <- as.matrix(dist_df)
  colnames(dist_matrix) <- NULL
  
  # annotation data.frame
  annotation_df <-
    sample_table_dds %>%
    dplyr::select(-sample_id)
  annotation_df$mutation <- factor(annotation_df$mutation, levels = c("delA", "P197A", "P216T",
                                                                      "A237V", "E176Q", "K70N",
                                                                      "G69V", "delHelicase"))
  
  # annotation colors
  ann_colors <- list(mutation = c("delA" = "#00B050", "P197A" = "#0000FF", "P216T" = "#FF0000", 
                                  "A237V" = "#FFC000", "E176Q" = "#00B0F0", "K70N" = "#A6A6A6", 
                                  "G69V" = "#C6E0B4", "delHelicase" = "#FF00FF"))
  
  # rownames annotation
  annotation_rownames <-
    rownames(dist_matrix) %>%
    str_remove_all(., "^s_|\\.PE$|\\.SE$") %>%
    str_replace_all(., "_", " ")
  
  # plot
  pheatmap::pheatmap(dist_matrix,
                     clustering_distance_rows = dist_df,
                     clustering_distance_cols = dist_df,
                     col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                     annotation_row = annotation_df,
                     annotation_colors = ann_colors, 
                     labels_row = annotation_rownames,
                     file = file.path(outpath, str_c("miRBase",
                                                     "plot", "dist_heatmap", "rlog", "mutation",
                                                     "labeled", "png", sep = ".")),
                     height = 10,
                     width = 14)
  
}

### FPM heatmap
if(TRUE){
  
  # annotation colors
  ann_colors <- list(mutation = c("delA" = "#00B050", "P197A" = "#0000FF", "P216T" = "#FF0000", 
                                  "A237V" = "#FFC000", "E176Q" = "#00B0F0", "K70N" = "#A6A6A6", 
                                  "G69V" = "#C6E0B4", "delHelicase" = "#FF00FF", 
                                  "null_mESCs" = "black"))
  
  # log transformed FPMs
  fpm_log_df <-
    fpm_tb %>%
    dplyr::select(gene_id, colnames(se_filt)) %>%
    dplyr::filter_at(.vars = vars(starts_with("s_")), .vars_predicate = any_vars(. > 0.5)) %>%
    dplyr::mutate_at(.vars = vars(starts_with("s_")), .funs = list(~ log2(. + 0.1))) %>%
    as.data.frame(.) %>%
    tibble::column_to_rownames(., var = "gene_id")
  
  # matrix for heatmap
  heatmap_matrix <- as.matrix(fpm_log_df)
  rownames(heatmap_matrix) <- NULL
  
  # annotation data.frame
  annotation_df <-
    sample_table_dds %>%
    dplyr::select(-c(sample_id))
  
  # sort rows and columns
  mat_cluster_cols <- sort_hclust(hclust(dist(t(heatmap_matrix))))
  mat_cluster_rows <- sort_hclust(hclust(dist(heatmap_matrix)))
  
  # plot
  pheatmap::pheatmap(heatmap_matrix,
                     col = viridis(10),
                     annotation_col = annotation_df,
                     annotation_colors = ann_colors, 
                     cluster_cols = mat_cluster_cols,
                     cluster_rows = mat_cluster_rows,
                     file = file.path(outpath, str_c("miRBase",
                                                     "plot", "FPM_heatmap", "log", "mutation",
                                                     "labeled", "png", sep = ".")),
                     height = 15,
                     width = 10)
  
}



### FPKM correlation
# FPKM values
fpm_corr_df <-
  fpm_tb %>%
  dplyr::select(colnames(se_filt))

### correlation matrix plot
# create plot
cor_pairs <- GGally::ggpairs(fpm_corr_df, diag = "blank")

# limit axis on all plots
for(i in 2:cor_pairs$nrow) {
  for(j in 1:(i - 1)) {
    cor_pairs[i, j] <-
      cor_pairs[i, j] +
      scale_x_continuous(limits = c(0, 500)) +
      scale_y_continuous(limits = c(0, 500))
  }
}

# add themes
cor_pairs <-
  cor_pairs +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

# save plot
png(filename = file.path(outpath, str_c("miRBase",
                                        "plot", "correlation_FPM", "log", "mutation",
                                        "labeled", "Dicer_KO", "png", sep = ".")), 
    width = 15, height = 15, units = "in", res = 300)
print(cor_pairs)
dev.off()
