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

# sort clusters
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

experiment='DicerX_mESC'
single_end=TRUE
threads=1
mapped_path='/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2018_Mar_Jun.mESC/Data/Mapped/STAR_mm10'
documentation_path='/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2018_Mar_Jun.mESC/Data/Documentation'
features_coordinates='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.gtf.gz'
features_name='ensembl.93.GRCm38.p6.20180919.UCSCseqnames'
genes_info_path='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv'
grouping_variables='genotype_clean'
results_groups='NULL'
protein_coding_only='no'
exploratory_analysis='yes'
interactive_plots='yes'
counts_path='./ensembl.93.GRCm38.p6.20180919.UCSCseqnames.counts.txt'
lfc_cut='2'
padj_cut='0.05'

# create and set outpath
outpath <- file.path(getwd(), "PCA_plots")
dir.create(outpath)

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# reads stats path
reads_stats_path <- list.files(mapped_path, ".*\\.stats_and_tracks\\.csv", full.names = T)

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.)))

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
### prepare tables
# get feature coordinates
features_tb <-
  counts_tb %>%
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>%
  as.data.table(.)

# get gene_id of protein coding genes
protein_genes <-
  genes_info %>%
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

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


### filter samples
# set list of chosen samples
sample_list <- list(DX_1_to_6 = c("s_ESC_DX_i1_JME1.SE", "s_ESC_DX_i2_JME2.SE", "s_ESC_DX_i3_JM7D1.SE", 
                                  "s_ESC_DX_i4_JM7D2.SE", "s_ESC_DX_i5_JM71.1.SE", "s_ESC_DX_i6_JM71.2.SE"), 
                    DXII_1_to_8 = c("s_ESC_DXII_i1.SE", "s_ESC_DXII_i2.SE", "s_ESC_DXII_i3.SE", 
                                    "s_ESC_DXII_i4.SE", "s_ESC_DXII_i5.SE", "s_ESC_DXII_i6.SE", 
                                    "s_ESC_DXII_i7.SE", "s_ESC_DXII_i8.SE"), 
                    DXII_3_to_8 = c("s_ESC_DXII_i3.SE", "s_ESC_DXII_i4.SE", "s_ESC_DXII_i5.SE", 
                                    "s_ESC_DXII_i6.SE", "s_ESC_DXII_i7.SE", "s_ESC_DXII_i8.SE"), 
                    DX_1_to_8 = c("s_ESC_DX_i1_JME1.SE", "s_ESC_DX_i2_JME2.SE", "s_ESC_DX_i3_JM7D1.SE", 
                                  "s_ESC_DX_i4_JM7D2.SE", "s_ESC_DX_i5_JM71.1.SE", "s_ESC_DX_i6_JM71.2.SE", 
                                  "s_ESC_DX_i7_R1BD1.SE", "s_ESC_DX_i8_R1BD2.SE"), 
                    all_samples = c("s_ESC_DX_i1_JME1.SE", "s_ESC_DX_i2_JME2.SE", "s_ESC_DX_i3_JM7D1.SE", 
                                    "s_ESC_DX_i4_JM7D2.SE", "s_ESC_DX_i5_JM71.1.SE", "s_ESC_DX_i6_JM71.2.SE", 
                                    "s_ESC_DX_i7_R1BD1.SE", "s_ESC_DX_i8_R1BD2.SE", "s_ESC_DXII_i1.SE", 
                                    "s_ESC_DXII_i2.SE", "s_ESC_DXII_i3.SE", "s_ESC_DXII_i4.SE", 
                                    "s_ESC_DXII_i5.SE", "s_ESC_DXII_i6.SE", "s_ESC_DXII_i7.SE", 
                                    "s_ESC_DXII_i8.SE"))

### loop through samples
purrr::map(names(sample_list), function(samples_name){
  
  # prepare sample table for DESeq colData
  sample_table_dds <-
    sample_table %>%
    dplyr::filter(sample_id %in% str_remove_all(colnames(se), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$")) %>% 
    dplyr::filter(sample_id != "s_ESC_DX_i3_JM7D1.SE") %>%
    dplyr::filter(sample_id %in% sample_list[[samples_name]]) %>% 
    dplyr::select(sample_id, clean_name, genotype_clean) %>% 
    dplyr::mutate(genotype_clean = factor(genotype_clean, levels = c("WT", "DicerX/X", "DicerX/X, low PKR", "DicerX/X + dPKR"))) %>% 
    as.data.frame(.) %>%
    set_rownames(., .$sample_id)
  
  # optionaly take only protein coding genes for analysis
  if(protein_coding_only == "yes"){
    
    # filter summarizedExperiment
    se_filt <- se[rownames(se) %in% protein_genes, ]
    
  }else{
    
    # don't filter
    se_filt <- se
    
  }
  
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
  dds <- DESeqDataSet(se_filt, design = ~genotype_clean)
  
  ### PCA plot
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
  
  
  ### plot
  # create bare plot object
  pca_plot <- ggplot(data = pca_tb, aes(x = PC1, y = PC2, label = sample_id))
  
  # color = first grouping variable
  pca_plot <-
    pca_plot +
    geom_point(aes_string(fill = "genotype_clean"), color = "black", size = 7.5, shape = 21) +
    scale_fill_manual(values = c("gray50", "red2", "#1a75ff", "blue3")) + 
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5)))
  
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
                                             "expl_plot", "PCA.PC1_PC2", "rlog", samples_name,
                                             "png", sep = ".")),
         plot = pca_plot, width = 12, height = 10)
  
  # add labels
  pca_plot <-
    pca_plot +
    geom_label_repel(aes(label = clean_name), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")
  
  # save labeled plot
  ggsave(filename = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                             "expl_plot", "PCA.PC1_PC2", "rlog", samples_name,
                                             "labeled", "png", sep = ".")),
         plot = pca_plot, width = 12, height = 10)
  
  # return
  return(samples_name)
  
})
