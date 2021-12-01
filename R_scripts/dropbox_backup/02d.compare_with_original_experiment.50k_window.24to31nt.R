### INFO:
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq/Analysis/expression.50k_window.24to31nt")

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
### in and out
# set inpath
inpath <- getwd()

# outpath
outpath <- getwd()

######################################################## READ DATA
### data from original experiment
# experiment 1 path
exp <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"
documentation_path <- file.path(exp, "Data/Documentation")
results_path <- file.path(exp, "Analysis/expression.50k_window.24to31nt")
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)
counts_path <- list.files(results_path, ".*\\.counts\\.txt$", full.names = T)
fpm_path <- list.files(results_path, ".*\\.FPM\\.csv", full.names = T)

# read counts from featureCounts
counts_tb_1 <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.))) %>%
  dplyr::mutate(Geneid = make.unique(Geneid))

# read FPM table
fpm_tb_1 <- readr::read_csv(fpm_path)

# read sample table
sample_table_1 <- data.table::fread(sample_table_path)


### data from re-sequencing experiment
# experiment 2 path
exp <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq"
documentation_path <- file.path(exp, "Data/Documentation")
results_path <- file.path(exp, "Analysis/expression.50k_window.24to31nt")
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)
counts_path <- list.files(results_path, ".*\\.counts\\.txt$", full.names = T)
fpm_path <- list.files(results_path, ".*\\.FPM\\.csv", full.names = T)
  
# read counts from featureCounts
counts_tb_2 <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.))) %>%
  dplyr::mutate(Geneid = make.unique(Geneid))

# read FPM table
fpm_tb_2 <- readr::read_csv(fpm_path)

# read sample table
sample_table_2 <- data.table::fread(sample_table_path)


### set grouping variables
# set grouping variables
grouping_variables <- c("genotype", "age")

######################################################## MAIN CODE
#### prepare data
## features
# get feature coordinates
features_tb <-
  counts_tb_1 %>%
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>%
  as.data.table(.)

## sample table
# clean sample table 1
sample_table_1 %<>% 
  dplyr::select(sample_id, stage, genotype, age) %>% 
  dplyr::filter(genotype != "Mov10l_HET")

# clean sample table 2
sample_table_2 %<>% 
  dplyr::select(sample_id, stage, genotype, age)

## count tables
# clean count table 1
counts_tb_1 %<>% 
  dplyr::select(-c(Chr:Length))

# clean count table 2
counts_tb_2 %<>% 
  dplyr::select(-c(Chr:Length))

# prepare sample table for DESeq colData
sample_table_dds <-
  rbind(sample_table_1, sample_table_2) %>%
  as.data.table(.) %>%
  .[sample_id %in% str_remove_all(c(colnames(counts_tb_1), colnames(counts_tb_2)), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$"), ] %>%
  .[, c("sample_id", grouping_variables), with = F] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)


## summarized experiment
# counts table
se <-
  dplyr::inner_join(counts_tb_1, counts_tb_2, by = "Geneid") %>%
  dplyr::rename(gene_id = Geneid) %>%
  dplyr::mutate_if(is.numeric, round, digits = 0) %>%
  as.data.frame(.) %>%
  set_rownames(., .$gene_id) %>%
  dplyr::select(-gene_id) %>%
  as.matrix(.)

# filter summarizedExperiment to include only chosen stage, set colData
se_filt <- se
se_filt <- SummarizedExperiment(list(counts = se_filt))
colnames(se_filt) <- str_remove_all(colnames(se_filt), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$")
se_filt <- se_filt[, colnames(se_filt)[match(rownames(sample_table_dds), colnames(se_filt))]]

# check if colnames of assay match rownames of sample table DataFrame
if(all(colnames(se_filt) == rownames(sample_table_dds))){
  
  # set sample table as colData
  colData(se_filt) <- DataFrame(sample_table_dds)
  
}else{
  
  # stop script with warrning
  stop("Columns in assay are not matching row of sample table. Please check your data annotation")
  
}


# ## FPM data for plots
# # data for plots = log transformed counts
# log_df <-
#   fpm_tb %>%
#   dplyr::select(-coordinates) %>%
#   dplyr::filter_at(.vars = vars(starts_with("s_")), .vars_predicate = any_vars(. > 0.5)) %>%
#   dplyr::mutate_at(.vars = vars(starts_with("s_")), .funs = list(~ log2(. + 0.1))) %>%
#   as.data.frame(.) %>%
#   tibble::column_to_rownames(., var = "gene_id") %>%
#   as.matrix(.)

## DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~grouped_variables)


#### exploratory analysis
### PCA plot
## calculate
# data for PCA = rlog transformed counts
rlog_df <-
  vst(dds, blind = T) %>%
  assay(.)

# calculates pca
pca <-
  rlog_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes table for ggplot
pca_tb <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))


## plot
# create bare plot object
pca_plot <- ggplot(data = pca_tb, aes(x = PC1, y = PC2, label = sample_id))

# if there is only one grouping variable use only color, if there is more use also a shape
if(length(grouping_variables) == 1){
  
  # color = first grouping variable
  pca_plot <-
    pca_plot +
    geom_point(aes_string(color = grouping_variables[1], fill = grouping_variables[1]), size = 5, shape = 21) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5)))
  
}else{
  
  # color = first grouping variable, shape = second grouping variable
  pca_plot <-
    pca_plot +
    geom_point(aes_string(color = grouping_variables[1], fill = grouping_variables[1], shape = grouping_variables[2]), size = 7.5) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5)),
           shape = guide_legend(override.aes = list(size = 5)))
  
}

# add labels, themes and save plot
pca_plot <-
  pca_plot +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# # save plot
# ggsave(filename = file.path(outpath, str_c("MesAur1.5k_windows",
#                                            "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
#                                            "compare", "png", sep = ".")),
#        plot = pca_plot, width = 12, height = 10)

# add labels
pca_plot <-
  pca_plot +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save labeled plot
ggsave(filename = file.path(outpath, str_c("MesAur1.50k_windows",
                                           "24to31nt",
                                           "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                           "labeled", "compare", "png", sep = ".")),
       plot = pca_plot, width = 12, height = 10)


### distance heatmap
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
  dplyr::select(-c(sample_id, grouped_variables))

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
                   labels_row = annotation_rownames,
                   file = file.path(outpath,
                                    str_c("MesAur1.50k_windows",
                                          "24to31nt",
                                          "plot", "dist_heatmap", "rlog", str_c(grouping_variables, collapse = "_"),
                                          "labeled", "compare", "png", sep = ".")),
                   height = 10,
                   width = 14)


# ### FPM correlation
# # FPM values
# fpm_corr_df <-
#   left_join(fpm_tb_1 %>% dplyr::select(-c(contains("HET") | contains("WT_21dpp"))), 
#             fpm_tb_2, by = c("gene_id", "coordinates")) %>%
#   dplyr::select(-c(gene_id, coordinates)) %>% 
#   .[, order(colnames(.))]
# 
# ### correlation matrix plot
# # create plot
# cor_pairs <- GGally::ggpairs(fpm_corr_df, diag = "blank")
# 
# # limit axis on all plots
# for(i in 2:cor_pairs$nrow) {
#   for(j in 1:(i - 1)) {
#     cor_pairs[i, j] <-
#       cor_pairs[i, j] +
#       scale_x_continuous(limits = c(0, 500)) +
#       scale_y_continuous(limits = c(0, 500))
#   }
# }
# 
# # add themes
# cor_pairs <-
#   cor_pairs +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # save plot
# png(filename = file.path(outpath,
#                          str_c("MesAur1.5k_windows",
#                                "plot", "correlation_FPM", "compare", "png", sep = ".")),
#     width = 20, height = 20, units = "in", res = 300)
# print(cor_pairs)
# dev.off()
