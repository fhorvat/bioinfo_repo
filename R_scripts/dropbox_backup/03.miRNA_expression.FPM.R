### INFO:
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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

######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

experiment='DicerX_embryos'
single_end=TRUE
threads=1
mapped_path='/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Data/Mapped/STAR_mm10'
documentation_path='/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Data/Documentation'
features_coordinates='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/miRBase.22.mm10.20181605.gff3'
features_name='miRBase.22.mm10.20181605'
genes_info_path=''
grouping_variables='genotype'
results_groups='DicerX_KO,DicerX_WT DicerX_HET,DicerX_WT DicerX_KO,DicerX_HET'
protein_coding_only='no'
exploratory_analysis='yes'
interactive_plots='yes'
counts_path='./miRBase.22.mm10.20181605.counts.txt'

# create and set outpath
outpath <- file.path(getwd(), str_c("results.", features_name))
dir.create(outpath)

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# reads stats path
reads_stats_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# FPM path
fpm_path <- list.files(inpath, str_c(features_name, "\\.FPM\\.csv$"), full.names = T)

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.))) %>% 
  dplyr::mutate(Geneid = make.unique(Geneid))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats table
reads_stats <-
  readr::read_delim(reads_stats_path, delim = "\t", col_names = c("sample_id", "library_size")) %>%
  dplyr::filter(!is.na(sample_id),
                str_detect(sample_id, "21to23nt$")) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.21to23nt$")) %>%
  as.data.table(.)

# read FPM table
fpm_tb <- readr::read_csv(fpm_path)

######################################################## MAIN CODE
### prepare tables
# get feature coordinates
features_tb <-
  counts_tb %>%
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>%
  as.data.table(.)

# prepare sample table for DESeq colData
sample_table_dds <-
  sample_table %>%
  as.data.table(.) %>%
  .[, c("sample_id", grouping_variables), with = F] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# data for plots = log transformed counts
log_df <-
  fpm_tb %>% 
  dplyr::select(-coordinates) %>% 
  dplyr::filter_at(.vars = vars(starts_with("s_")), .vars_predicate = any_vars(. > 0.5)) %>% 
  dplyr::mutate_at(.vars = vars(starts_with("s_")), .funs = list(~ log2(. + 0.1))) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_id") %>% 
  as.matrix(.)


### PCA plot ####
# calculates pca
pca <-
  log_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes table for ggplot
pca_tb <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         sample_id = colnames(log_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))

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

# save plot
ggsave(filename = file.path(outpath, str_c("miRBase",
                                           "plot", "PCA.PC1_PC2", "log", str_c(grouping_variables, collapse = "_"),
                                           "png", sep = ".")),
       plot = pca_plot, width = 12, height = 10)

# add labels
pca_plot <-
  pca_plot +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save labeled plot
ggsave(filename = file.path(outpath, str_c("miRBase",
                                           "plot", "PCA.PC1_PC2", "log", str_c(grouping_variables, collapse = "_"),
                                           "labeled", "png", sep = ".")),
       plot = pca_plot, width = 12, height = 10)



### FPM heatmap ####
# make matrix
heatmap_matrix <- as.matrix(log_df)
rownames(heatmap_matrix) <- NULL

# annotation data.frame
annotation_df <-
  sample_table_dds %>% 
  dplyr::select(-c(sample_id, grouped_variables))

# sort clusters
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

# sort rows and columns 
mat_cluster_cols <- sort_hclust(hclust(dist(t(heatmap_matrix))))
mat_cluster_rows <- sort_hclust(hclust(dist(heatmap_matrix)))

# plot
pheatmap::pheatmap(heatmap_matrix,
                   col = viridis(10),
                   annotation_col = annotation_df,
                   cluster_cols = mat_cluster_cols,
                   cluster_rows = mat_cluster_rows,
                   file = file.path(outpath,
                                    str_c("miRBase",
                                          "plot", "FPM_heatmap", "log", str_c(grouping_variables, collapse = "_"),
                                          "png", sep = ".")),
                   height = 15,
                   width = 10)



## MA plot
# data for plot
plot_df <- 
  fpm_tb %>% 
  dplyr::select(-coordinates) %>% 
  tidyr::pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "fpm") %>% 
  dplyr::left_join(., sample_table %>% dplyr::select(sample_id, genotype), by = "sample_id") %>% 
  dplyr::group_by(gene_id, genotype) %>% 
  dplyr::summarise(fpm_avg = mean(fpm)) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::pivot_wider(id_cols = gene_id, names_from = genotype, values_from = fpm_avg) %>% 
  dplyr::mutate(mean = ((WT + KO) / 2),
                lfc = log2(KO + 0.1) - log2(WT + 0.1), 
                arm = str_extract(gene_id, "3p|5p")) %>% 
  dplyr::filter(!is.na(arm)) %>% 
  dplyr::select(gene_id, mean, lfc, arm) 

# result limits
results_limits <-
  plot_df %>% 
  dplyr::summarise(x_limit = mean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = lfc %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# plot
ma_plot <-
  ggplot() +
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = arm), size = 5, shape = 20) +
  scale_x_log10(limits = c(0.001, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                     breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
  xlab("mean expression") +
  ylab(str_c("log2 fold change: ", "DicerX KO", " / ", "DicerX WT", "\n")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) +
  theme(legend.title = element_blank())

# # turns off axis titles and legend
# ma_plot <-
#   ma_plot +
#   theme(legend.position = "none") +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())

# save plot
ggsave(filename = file.path(outpath,
                            str_c("miRBase",
                                  "plot", "MA", "FPKM", str_c(grouping_variables, collapse = "_"),
                                  "png", sep = ".")),
       plot = ma_plot, width = 12, height = 10)


