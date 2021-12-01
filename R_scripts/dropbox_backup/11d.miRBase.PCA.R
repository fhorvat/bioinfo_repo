### INFO: counts reads over mature miRNAs (from miRbase .gff)
### DATE: Wed Nov 28 21:54:06 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(DESeq2)
library(rtracklayer)
library(xlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get PCA and variance of each PC
df2PCA <- function(df, sample_table, plot_name){

  # calculate PCA
  pca <-
    df %>%
    t(.) %>%
    stats::prcomp(.)
  
  # gets percent of variance for each principal component
  percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)
  
  # makes data.frame for ggplot
  plot_df <-
    tibble(PC1 = pca$x[, 1],
           PC2 = pca$x[, 2],
           sample_id = colnames(df)) %>%
    dplyr::left_join(sample_table , by = "sample_id") %>%
    dplyr::mutate(sample_id = str_remove_all(sample_id, "^s_|r") %>% str_replace_all(., "_", " "))
  
  # plot 
  pca_plot <- 
    ggplot(data = plot_df, aes(x = PC1, y = PC2, label = sample_id, fill = genotype, shape = transfection)) +
    geom_point(aes(fill = genotype), color = "black", size = 7.5) +
    scale_shape_manual(values  = c(21, 22)) +
    xlab(str_c("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(str_c("PC2: ", round(percentVar[2] * 100), "% variance")) +
    # guides(fill = FALSE, shape = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 23, size = 5)),
           shape = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
          axis.title.y = element_text(size = 15, vjust = - 0.2),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, str_c("PCAplot.miRBase.no_labels.", plot_name, ".png")), width = 12, height = 10)
  
  # add labels and save
  pca_plot <- 
    pca_plot + 
    ggrepel::geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", 
                              box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
    ggsave(filename = file.path(outpath, str_c("PCAplot.miRBase.with_labels.", plot_name, ".png")), width = 12, height = 10)
  
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# summarizedOverlaps .RDS path
se_path <- file.path(inpath, "miRbase.Eliska_mESC_MosIR.se.RDS")

# miRNA RPM path
mirna_rpm_path <- file.path(inpath, "mature_miRNA.Eliska_mESC_MosIR.expression.xlsx")

# sample table path
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Eliska_mESC_MosIR/Data/Documentation/Eliska_mESC_MosIR.sampleTable.csv"

######################################################## READ DATA
# read summarizeOverlaps .RDS
se <- readRDS(file = se_path)

# read library size
mirna_rpm <- xlsx::read.xlsx(file = mirna_rpm_path, sheetName = "allowed_0_mismatches.RPMs")

# read sample table
sample_table <- readr::read_csv(sample_table_path)

######################################################## MAIN CODE
### summarizedExperiment
# change column names of SE
colnames(se) <- str_remove(colnames(se), "\\.SE\\.mis_0\\.18to30nt\\.bam")

# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE$"), 
                group = str_c(genotype, "_", transfection)) %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# rearrange sample table so it corrresponds to colnames of SE
sample_table_dds <- sample_table_dds[match(colnames(se), rownames(sample_table_dds)), ]

# add colData to SE
colData(se) <- DataFrame(sample_table_dds)

# make DESeqDataSet
dds <- DESeqDataSet(se, design = ~group)

# data for PCA = rlog transformed counts
pca_plot <-
  rlog(dds, blind = T) %>%
  assay(.) %>% 
  df2PCA(df = ., sample_table = sample_table_dds, plot_name = "rlog") 

### RPM
pca_plot <- 
  mirna_rpm %>% 
  dplyr::select(-mirna_id, -coordinates) %>% 
  dplyr::mutate_all(funs(log(. + 1))) %>% 
  df2PCA(df = ., sample_table = sample_table_dds, plot_name = "RPM") 

