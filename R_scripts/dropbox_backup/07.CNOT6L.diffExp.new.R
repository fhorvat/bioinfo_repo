### INFO: CNOT6L differential expression analysis
### DATE: Sun Mar 11 04:54:44 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/review/diffExp")

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
library(plotly)
library(VennDiagram)
library(geneplotter)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/2018_paper/expression_analysis"

# set main outpath
outpath <- getwd()

# ENSEMBL annotated genes info path
ensembl_genes_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.89.GRCm38.p5.20180615.UCSCseqnames.geneInfo.csv"

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
    dplyr::mutate(comparison = str_c(filt_stage, ".KO_vs_WT")) %>%
    write_csv(., path = file.path(outpath, str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT.padj005", ".GRCm38.89.all.csv")))

  # write only significant results, padj < 0.05
  results_df %>%
    dplyr::filter(padj < 0.05) %T>%
    write_csv(., path = file.path(outpath, str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT.padj005", ".GRCm38.89.sign.csv")))
  
  # read results
  results_df <- readr::read_csv(file = file.path(outpath, str_c("diffExp.CNOT6L.", filt_stage, ".KO_vs_WT.padj005", ".GRCm38.89.all.csv")))
  
  ### DESeq2 MA plot
  # data for plot
  plot_df <- 
    results_df %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  sign = ifelse(padj < 0.05, "yes", "no"),
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, sign == "no", "not_sign"),
                  regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>%
    dplyr::arrange(regulation)
  
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
    ggtitle(str_c(filt_stage, " KO vs. WT")) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT.padj005", ".GRCm38.89.png")),
         plot = ma_plot, width = 10, height = 10)
  
  
  # ### FPKM MA plot
  # # data.frame for plot
  # plot_df <-
  #   results_df %>%
  #   dplyr::select(matches(filt_stage), padj, gene_id) %>%
  #   dplyr::rename_at(.vars = vars(matches("KO|WT")), .funs = funs(c("lfc_KO", "lfc_WT"))) %>%
  #   dplyr::mutate(mean = dplyr::select(., lfc_KO, lfc_WT) %>% rowMeans(., na.rm = T),
  #                 lfc = (log2(lfc_KO + 1) - log2(lfc_WT + 1)),
  #                 padj = replace(padj, is.na(padj), 1),
  #                 sign = ifelse(padj < 0.05, "yes", "no"),
  #                 regulation = ifelse(lfc > 0, "up", "down"),
  #                 regulation = replace(regulation, sign == "no", "not_sign"),
  #                 regulation = factor(regulation, levels = c("not_sign", "up", "down"))) %>%
  #   dplyr::arrange(regulation)
  # 
  # # plot
  # ma_plot <-
  #   ggplot(data = plot_df, aes(x = mean, y = lfc, color = regulation)) +
  #   geom_point(size = 3, shape = 20) +
  #   scale_x_log10(limits = c(1e-03, 1e5),
  #                 breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  #   scale_y_continuous(limits = c(-7, 5),
  #                      breaks = c(-7:5)) +
  #   scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
  #   guides(color = FALSE) +
  #   ggtitle(str_c(filt_stage, " KO vs. WT")) +
  #   theme_bw() +
  #   theme(axis.title.x = element_blank(),
  #         axis.title.y = element_blank()) +
  #   theme(axis.text.x = element_text(size = 15),
  #         axis.text.y = element_text(size = 15),
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank())
  # 
  # # save plot
  # ggsave(filename = file.path(outpath, str_c("MAplot.CNOT6L.", filt_stage, ".KO_vs_WT.padj005", ".FPKM.GRCm38.89.png")),
  #        plot = ma_plot, width = 10, height = 10)
  
}

### Venn diagrams
# read significantly up- and down-regulated genes
signif_genes <- list.files(outpath, pattern = "*.signif.csv")
signif_genes <- 
  lapply(signif_genes, readr::read_csv) %>% 
  magrittr::set_names(str_extract(signif_genes, "GV|MII|1C"))

### create list of gene_id
# down-regulated
down_genes <- lapply(signif_genes, function(x){x %>% 
    dplyr::filter(log2FoldChange < 0) %$% gene_id})

# up-regulated
up_genes <- lapply(signif_genes, function(x){x %>% 
    dplyr::filter(log2FoldChange > 0) %$% gene_id})

### plot
# down-regulated
png(file = file.path(outpath, "Venn.CNOT6L.significant.downregulated.KO_vs_WT.padj005.GRCm38.89.png"), 
    width = 1000, height = 1000, units = "px", type = "cairo")
down_plot <- venn.diagram(x = down_genes, 
                          filename = NULL, 
                          fill = c("blue", "red", "green"), 
                          alpha = c(0.5, 0.5, 0.5), 
                          cex = 4,
                          cat.cex = 2,
                          category.names = names(down_genes), 
                          main = "Downregulated", 
                          main.cex = 2)
grid.draw(down_plot)
dev.off()

# up-regulated 
png(file = file.path(outpath, "Venn.CNOT6L.significant.upregulated.KO_vs_WT.padj005.GRCm38.89.png"), 
    width = 1000, height = 1000, units = "px", type = "cairo")
up_plot <- venn.diagram(x = up_genes, 
                        filename = NULL, 
                        fill = c("blue", "red", "green"), 
                        alpha = c(0.5, 0.5, 0.5), 
                        cex = 4,
                        cat.cex = 2,
                        category.names = names(up_genes), 
                        main = "Upregulated", 
                        main.cex = 2)
grid.draw(up_plot)
dev.off()










