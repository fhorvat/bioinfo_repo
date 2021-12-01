### INFO: CNOT6L differential expression analysis
### DATE: Sun Mar 11 04:54:44 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

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

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.geneInfo.csv"

######################################################## READ DATA
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

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
se_filt <- se[rownames(se) %in% protein_genes, ]
colData(se_filt) <- DataFrame(sample_table_dds)
se_filt$genotype <- factor(se_filt$genotype, levels = c("WT", "KO"))
se_filt$group <- factor(se_filt$group, levels = c("GV_WT", "GV_KO", "MII_WT", "MII_KO", "1C_WT", "1C_KO"))

# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~group)

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
    tibble::rownames_to_column(., var = "gene_id") %>%
    as.tibble(.) %>%
    dplyr::arrange(padj) %>%
    dplyr::left_join(fpkm_df %>% dplyr::select(gene_id, matches(filt_stage)), by = "gene_id") %>%
    dplyr::left_join(ensembl_genes_info_filt, by = "gene_id") %>%
    dplyr::mutate(comparison = str_c(filt_stage, ".KO_vs_WT")) %T>%
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv")))
  
  # write only significant results, padj < 0.1
  results_df %>%
    dplyr::filter(padj < 0.1) %>%
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.signif.csv")))
  
  ### MA plot
  # read results 
  results_df <- read_csv(file = file.path(outpath, "results", str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.all.csv")))
  
  # data for plot
  plot_df <- 
    results_df %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
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
  ggsave(filename = file.path(outpath, "results", str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT", ".GRCm38.89.png")),
         plot = ma_plot, width = 10, height = 10)
  
}
