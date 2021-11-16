### INFO: 
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY

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

######################################################## PATH VARIABLES
### experiment
# set experiment name
experiment <- "Lnc1_KO"
experiment_name <- "Lnc1_KO"

### working dir
# set working directory 
setwd(file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/Analysis/lnc1_paper/diffExp/lnc1_KO"))


### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


### other experiment paths
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec")

# mapped path
mapped_path <- file.path(base_path, "Data/Mapped/STAR_mm10")

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- inpath


### documentation
# set ensembl version
ensembl_version <- 93

# sample table path
sample_table_path <- file.path(documentation_path, "lnc1_KO.RNAseq.20181211.sampleTable.clean.csv")

# stats and tracks path
stats_and_tracks_path <- list.files(path = mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.se\\.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.FPKM_mean\\.csv$"), full.names = T)


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

# read stats and tracks table
stats_and_tracks <- data.table::fread(stats_and_tracks_path)

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
# join sample table with stats and tracks
sample_table[stats_and_tracks, on = "sample_id", `:=`(library_size = genome.mapped_minus_rDNA)]

# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter ensembl genes info
genes_info_tidy <- 
  genes_info %>% 
  dplyr::mutate(cooridnates = str_c(seqnames, ":", start, "-", end, "|", strand)) %>% 
  dplyr::select(-c(seqnames:strand))

# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  tidyr::unite(col = "group", genotype, stage, sep = "_", remove = F) %>% 
  dplyr::select(sample_id, genotype, stage, group) %>% 
  as.data.frame(.) %>% 
  set_rownames(., .$sample_id) 


### DESeq2
# set stages
stages <- c("GV", "MII")

# create results workbooks
wb_all <- openxlsx::createWorkbook()
wb_significant <- openxlsx::createWorkbook()


### get all results
purrr::map(stages, function(stage){
  
  
  ### prepare sample table and summarizedExperiment
  # filter sample table
  sample_table_filt <- sample_table_dds[sample_table_dds$stage == stage, ]
  
  # filter summarizedExperiment to include only chosen stage and protein coding genes, set colData
  se_filt <- se
  se_filt <- se_filt[rownames(se_filt) %in% c(protein_genes, "ENSMUSG00000110001"), ] 
  colnames(se_filt) <- str_remove(colnames(se_filt), "\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$")
  se_filt <- se_filt[, colnames(se_filt)[match(rownames(sample_table_filt), colnames(se_filt))]]
  
  # check if colnames of assay match rownames of sample table DataFrame
  if(all(colnames(se_filt) == rownames(sample_table_filt))){
    
    # set sample table as colData
    colData(se_filt) <- DataFrame(sample_table_filt)
    
  }else{
    
    # stop script with warrning
    stop("Columns in assay are not matching row of sample table. Please check your data annotation")
    
  }
  
  
  ### do DESeq and results
  # make DESeqDataSet
  dds <- 
    DESeqDataSet(se_filt, design = ~genotype) %>% 
    DESeq(.)
  
  # shrink results
  dds_shrink <- lfcShrink(dds, contrast = c("genotype", "Null", "WT"))
  
  # get results table
  results_df <-
    dds_shrink %>%
    as_tibble(., rownames = "gene_id") %>%
    dplyr::arrange(padj) %>%
    dplyr::left_join(fpkm_tb %>% dplyr::select_at(.vars = vars(gene_id, contains(stage), gene_name:gene_description)), 
                     by = "gene_id") %>%
    dplyr::mutate(comparison = str_c(stage, ".lnc1_Null_vs_WT")) %>% 
    setnames(x = ., 
             old = colnames(.)[str_detect(colnames(.), stage)], 
             new = str_c(colnames(.)[str_detect(colnames(.), stage)], ".FPKM"))

  # # write results to sheet
  # openxlsx::addWorksheet(wb = wb_all, sheetName = str_c(stage, ".lnc1_Null_vs_WT"))
  # openxlsx::writeData(wb = wb_all,
  #                     sheet = str_c(stage, ".lnc1_Null_vs_WT"),
  #                     x = results_df)

  
  ### write only significant results, padj < 0.1
  # filter table
  results_df_sign <- 
    results_df %>% 
    dplyr::filter(padj < 0.1)
  
  # # check and write
  # if(nrow(results_df_sign) > 0){
  # 
  #   # write data
  #   openxlsx::addWorksheet(wb = wb_significant, sheetName = str_c(stage, ".lnc1_Null_vs_WT"))
  #   openxlsx::writeData(wb = wb_significant,
  #                       sheet = str_c(stage, ".lnc1_Null_vs_WT"),
  #                       x = results_df_sign)
  # 
  # }
  
  ### MA plot - ggplot2
  # data for plot
  plot_df <- 
    results_df %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, coordinates) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, padj > 0.1, "no")) %>%
    dplyr::mutate(regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
    dplyr::mutate(mito = ifelse(str_detect(coordinates, "chrM"), "mitochondrion", "other"),
                  mito = factor(mito, levels = c("other", "mitochondrion"))) %>%
    dplyr::arrange(regulation, mito)
  
  # get annotation data
  annotation_df <- 
    plot_df %>% 
    dplyr::filter(gene_id %in% c("ENSMUSG00000055839", "ENSMUSG00000110001", "ENSMUSG00000057534"))
  
  # plot
  ma_plot <- 
    ggplot() + 
    geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation, shape = mito), size = 5) +
    geom_point(data = annotation_df, aes(x = mean, y = lfc, shape = mito), color = "black", alpha = 1, size = 5) +
    scale_x_log10(limits = c(0.1, 1e6), 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_continuous(limits = c(-3, 3), breaks = c(-3:3)) +
    scale_colour_manual(values = c(no = "gray50", down = "#1a75ff", up = "red3")) +
    scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
    scale_shape_manual(values = c(other = 20, mitochondrion = 17)) +
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
  ggsave(filename = file.path(outpath, str_c("diffExp", 
                                             "ensembl", ensembl_version, 
                                             stage, "lnc1_Null_vs_WT", 
                                             "DESeq2.protein_coding.MA_plot.png", sep = ".")),
         plot = ma_plot, width = 10, height = 10)
  
  # add lables to 3 genes
  ma_plot <- 
    ma_plot + 
    geom_label_repel(data = annotation_df, aes(x = mean, y = lfc, label = gene_name), 
                     fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("diffExp", 
                                             "ensembl", ensembl_version, 
                                             stage, "lnc1_Null_vs_WT", 
                                             "DESeq2.protein_coding.MA_plot.labels.png", sep = ".")),
         plot = ma_plot, width = 10, height = 10)
  
})

# save workbooks to outdir
openxlsx::saveWorkbook(wb = wb_all, 
                       file = file.path(outpath, str_c("diffExp", 
                                                       "ensembl", ensembl_version,
                                                       "lnc1_Null_vs_WT", 
                                                       "DESeq2.protein_coding.all_results.xlsx", sep = ".")), 
                       overwrite = TRUE)

openxlsx::saveWorkbook(wb = wb_significant, 
                       file = file.path(outpath, str_c("diffExp", 
                                                       "ensembl", ensembl_version,
                                                       "lnc1_Null_vs_WT",
                                                       "DESeq2.protein_coding.significant_results.xlsx", sep = ".")), 
                       overwrite = TRUE)


