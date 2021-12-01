### INFO: 
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/diffExp/Stein_2015_PLoSGenet_GSE57514")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### experiment
# set experiment name
experiment <- "Stein_2015_PLoSGenet_GSE57514"

# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10_new")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis")


### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.se.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.*.avgFPKM.csv$"), full.names = T)


### genome
# genome path
genome_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# summarizedExperiment 
se <- readRDS(se_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
### prepare data
# filter ensembl genes info
genes_info_tidy <- 
  genes_info %>% 
  dplyr::mutate(cooridnates = str_c(seqnames, ":", start, "-", end, "|", strand)) %>% 
  dplyr::select(-c(seqnames:strand))

# # get gene_id of protein coding genes
# protein_genes <- 
#   genes_info_tidy %>% 
#   dplyr::filter(gene_biotype == "protein_coding") %$%
#   gene_id


### DESeq2
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size) %>% 
  dplyr::filter(str_detect(genotype, "Dicer")) %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
# se_filt <- se[rownames(se) %in% protein_genes, ]
se_filt <- se
colnames(se_filt) <- str_remove(colnames(se_filt), ".genome.Aligned.sortedByCoord.out.bam|.total.bam")
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
dds <- DESeqDataSet(se_filt, design = ~genotype)

dds_norm <- 
  dds %>% 
  estimateSizeFactors(.) %>% 
  counts(., normalized = TRUE)


# run DESeq
dds_deseq <- DESeq(dds)

# get results, shrink logFC
dds_shrink <- lfcShrink(dds_deseq, contrast = c("genotype", "Dicer_KO", "Dicer_WT"))

# get results table
results_df <-
  dds_shrink %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  as.tibble(.) %>%
  dplyr::arrange(padj) %>%
  dplyr::left_join(fpkm_tb %>% dplyr::select(gene_id, contains("Dicer")), by = "gene_id") %>%
  dplyr::left_join(genes_info_tidy, by = "gene_id") %>%
  dplyr::mutate(comparison = "Dicer.KO_vs_WT")
  # write_csv(., path = file.path(outpath, str_c("diffExp", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "all_biotype.all.csv", sep = ".")))

# write only significant results, padj < 0.1
results_df %>%
  dplyr::filter(padj < 0.1)
  # write_csv(., path = file.path(outpath, str_c("diffExp", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "all_biotype.signif.csv", sep = ".")))


### DESeq2 MA plot
# data for plot
plot_df <- 
  results_df %>% 
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, padj > 0.1, "no")) %>%
  dplyr::mutate(regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
  dplyr::arrange(regulation)

# get annotation data
annotation_df <- 
  plot_df %>% 
  dplyr::filter(gene_id %in% c("ENSMUSG00000055839", "ENSMUSG00000110001", "ENSMUSG00000057534")) %T>%
  readr::write_csv(., "Sirena_Elob_Elobl.diffExp.Dicer.KO_vs_WT.Stein.csv")

# plot
ma_plot <- 
  ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 2.5, shape = 20) +
  geom_point(data = annotation_df, aes(x = mean, y = lfc), shape = 20, colour = "black", fill = "black", size = 2.5) +
  scale_x_log10(limits = c(0.1, 1e7), 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  scale_y_continuous(limits = c(-6, 6), breaks = c(-6:6)) +
  scale_colour_manual(values = c(no = "gray50", down = "#1a75ff", up = "red2")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "none")

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "all_biotype.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)

# add lables to 3 genes
ma_plot <- 
  ma_plot + 
  geom_label_repel(data = annotation_df, aes(x = mean, y = lfc, label = gene_name), 
                   fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save plot
ggsave(filename = file.path(outpath, str_c("MAplot", experiment, "Dicer_KO_vs_WT", "ensembl", ensembl_version, "all_biotype.labels.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)