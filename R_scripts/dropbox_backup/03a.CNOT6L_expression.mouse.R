### INFO: get expression of all genes in CNOT6L data
### DATE: Sat Mar 10 19:01:11 2018
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

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# reduced exons path
exons_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.reducedExons.RDS"

# stats and tracks, bam paths
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10_noMultimapFilter"

# info about chosen genes path
genes_info_path <- file.path(outpath, "CNOT.GRCm38.89.20180305.geneInfo.csv")

# aditional info about genes path
ensembl_genes_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.geneInfo.csv"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read ENSEMBL reduced exons 
exons_gr <- readRDS(file = exons_path)

# read stats and tracks table
stats_df <- readr::read_csv(file = list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T))

# read info about chosen genes
genes_info <- readr::read_csv(genes_info_path)

# read info about all ensembl annotated genes
ensembl_genes_info <- readr::read_csv(ensembl_genes_path)

######################################################## MAIN CODE
# get total length of all exons for each transcript
exons_width <- 
  width(exons_gr) %>%
  sum(.) %>% 
  tibble(gene_id = names(.), width = .)

# construct sample table
sample_table <- 
  tibble(bam_path = list.files(path = mapped_path, pattern = "*.bam$", full.names = T)) %>% 
  dplyr::mutate(sample_id = str_replace(basename(bam_path), ".Aligned.sortedByCoord.out.bam", "")) %>% 
  dplyr::left_join(stats_df %>% dplyr::select(sample_id, library_size = mapped_genome_except_rDNA), by = "sample_id") %>% 
  dplyr::mutate(stage = str_extract(sample_id, "1C|MII|GV"), 
                genotype = str_extract(sample_id, "WT|KO"))

# # get count of reads, save summarizedExperiment as RDS
# bamfiles <- Rsamtools::BamFileList(sample_table$bam_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = exons_gr,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"))

# read summarizedExperiment from RDS file
se <- readRDS(file = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.summarizedOverlaps.RDS"))

# ### FPKM
# # get data.frame of counts, transform to FPKM
# fpkm_df_avg <-
#   assay(se) %>%
#   as.data.frame(.) %>%
#   tibble::rownames_to_column(., var = "gene_id") %>%
#   as.tibble(.) %>%
#   set_colnames(., str_replace(colnames(.), ".Aligned.sortedByCoord.out.bam", "")) %>%
#   tidyr::gather(key = sample_id, value = counts, -gene_id) %>%
#   dplyr::left_join(., sample_table %>% dplyr::select(sample_id, library_size), by = "sample_id") %>%
#   dplyr::left_join(., exons_width, by = "gene_id") %>%
#   dplyr::mutate(library_size = round(library_size / 1E6, 6),
#                 width = round(width / 1E3, 3),
#                 fpm = (counts / library_size),
#                 fpkm = (fpm / width)) %>%
#   dplyr::select(gene_id, sample_id, fpkm) %>%
#   dplyr::mutate(stage = str_replace_all(sample_id, "_r[1-3]{1}", "")) %>%
#   dplyr::group_by(gene_id, stage) %>%
#   dplyr::summarise(average_fpkm = mean(fpkm)) %>%
#   dplyr::ungroup(.) %>%
#   tidyr::spread(key = stage, value = average_fpkm) %T>%
#   readr::write_csv(., path = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

# fpkm_df %>% 
#   dplyr::select(gene_id, contains("MII|1C")) %>%
#   dplyr::filter(gene_id == "ENSMUSG00000034724") %>% 
#   as.data.frame(.)
# 
# fpkm_df_avg %>% 
#   dplyr::select(gene_id, matches("s_MII.*|s_1C.*")) %>%
#   dplyr::filter(gene_id == "ENSMUSG00000034724") %>% 
#   as.data.frame(.) %>% 
#   dplyr::mutate(log2FC_1C = log2(s_1C_KO) - log2(s_1C_WT), 
#                 log2FC_MII = log2(s_MII_KO) - log2(s_MII_WT))

# read FPKM table
fpkm_df <- readr::read_csv(file = file.path(outpath, "ensembl.GRCm38.89.CNOT6L.avgFPKM.csv"))

# # get FPKM of chosen genes
# fpkm_df %>%
#   dplyr::right_join(., genes_info %>% dplyr::select(gene_id, gene_name), by = "gene_id") %>%
#   dplyr::select(gene_name, everything(), -gene_id) %>%
#   dplyr::arrange(gene_name) %T>%
#   readr::write_csv(., path = file.path(outpath, "CNOT.GRCm38.89.CNOT6L.avgFPKM.csv"))

### differential expression analysis
# filter sample table
sample_table_filt <- 
  sample_table %>% 
  dplyr::select(-bam_path, -library_size)

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


# go through stages
for(filt_stage in unique(sample_table_filt$stage)){
  
  # filter FPKM table
  fpkm_df_filt <-
    fpkm_df %>%
    dplyr::select(gene_id, contains(filt_stage)) %>%
    data.table::setnames(., -1, str_c("mean_FPKM_", colnames(.)[2:ncol(.)]))
  
  # filter sample table
  sample_table_stage <-
    sample_table_filt %>%
    dplyr::filter(stage == filt_stage) %>%
    as.data.frame(.) %>%
    set_rownames(., .$sample_id)
  
  # filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
  se_filt <- se[, str_detect(colnames(se), filt_stage)]
  se_filt <- se_filt[rownames(se_filt) %in% protein_genes, ]
  colData(se_filt) <- DataFrame(sample_table_stage)
  se_filt$genotype <- factor(se_filt$genotype, levels = c("WT", "KO"))
  
  # run DESeq
  dds <-
    DESeqDataSet(se_filt, design = ~genotype) %>%
    DESeq(.)
  
  # get results, shrink logFC
  dds_shrink <- lfcShrink(dds, coef = 2)
  saveRDS(dds_shrink, file = str_c("diffExp.CNOT6L.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.lfcShrink.RDS"))

  # read lfcShrink transformed results
  dds_shrink <- readRDS(file = str_c("diffExp.CNOT6L.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.lfcShrink.RDS"))
  
  # get results table
  results_df <-
    dds_shrink %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(., var = "gene_id") %>%
    as.tibble(.) %>%
    dplyr::arrange(padj) %>%
    dplyr::left_join(fpkm_df_filt, by = "gene_id") %>%
    dplyr::left_join(ensembl_genes_info_filt, by = "gene_id") %T>%
    write_csv(., path = file.path(outpath, "results",
                                  str_c("diffExp.CNOT6L.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.all.csv")))
  
  # write only significant results, padj < 0.1
  results_df %>%
    dplyr::filter(padj < 0.1) %>%
    write_csv(., path = file.path(outpath, "results",
                                  str_c("diffExp.CNOT6L.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.signif.csv")))
  
  # read results
  results_df <- read_csv(file = file.path(outpath, "results", 
                                          str_c("diffExp.CNOT6L.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.all.csv")))
  
  ### MA plot
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
    scale_y_continuous(limits = c(-5.1, 2.5), 
                       breaks = c(-5:2), 
                       labels = c("", "-4", "", "-2", "", "", "", "2")) +
    scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
    guides(alpha = FALSE, color = FALSE) +
    # xlab("") +
    # ylab("") +
    ggtitle(str_c(filt_stage, " ", str_replace_all(resultsNames(dds)[2], "_", " ")))
    theme_bw() +
    # theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
    #       axis.title.y = element_text(size = 15, vjust = 0.3)) +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", 
                              str_c("MAplot.CNOT6L.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.noShrink.png")),
         plot = ma_plot, width = 10, height = 10)
  
}
