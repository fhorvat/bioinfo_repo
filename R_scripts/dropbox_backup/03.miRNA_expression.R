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
reads_stats_path <- file.path(mapped_path, "library_sizes.txt")

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

# read stats and tracks table
reads_stats <- data.table::fread(reads_stats_path)
data.table::setnames(reads_stats, c("sample_id", "library_size"))
reads_stats[, sample_id := str_remove(sample_id, "\\.21to23nt")]

# read FPM table
fpm_tb <- readr::read_csv(fpm_path)

######################################################## MAIN CODE
### prepare tables
# get feature coordinates
features_tb <-
  counts_tb %>%
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>%
  as.data.table(.)

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

# join sample table with stats and tracks
sample_table[reads_stats, on = "sample_id", `:=`(library_size = library_size)]

# create vector of plotly symbols in ggplot shape order
ploty_symbols <- c("square", "circle", "cross", "x", "diamond", "square-open", "circle-open", "diamond-open")


### filter samples
# prepare sample table for DESeq colData
sample_table_dds <-
  sample_table %>%
  as.data.table(.) %>%
  .[, c("sample_id", grouping_variables), with = F] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# don't filter protein coding genes since this is miRNAs
se_filt <- se

# filter summarizedExperiment to include only chosen stage, set colData
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


### DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~grouped_variables)


####### EXPLORATORY ANALYSIS
if(exploratory_analysis == "yes"){
  
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
  
  # makes table for ggplot
  pca_tb <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           sample_id = colnames(rlog_df)) %>%
    dplyr::left_join(sample_table_dds , by = "sample_id") %>%
    dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))
  
  
  ### plot
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
  ggsave(filename = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                             "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                             "png", sep = ".")),
         plot = pca_plot, width = 12, height = 10)
  
  # add labels
  pca_plot <-
    pca_plot +
    geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")
  
  # save labeled plot
  ggsave(filename = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                             "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                             "labeled", "png", sep = ".")),
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
                                      str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                            "plot", "dist_heatmap", "rlog", str_c(grouping_variables, collapse = "_"),
                                            "png", sep = ".")),
                     height = 10,
                     width = 14)
  
}

