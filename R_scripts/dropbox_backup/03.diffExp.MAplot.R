### INFO: Do the differential expression analysis for mESC and oocytes sequenced in February 2018
### DATE: Mon Mar 05 16:05:23 2018
### AUTHOR: Filip Horvat

rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Analysis/expression")

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
library(plotly)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "GffToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Documentation/mESC_oocytes_2018.sample_table.csv"

# genes info
genes_info_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Mus_musculus.GRCm38.89.20180305.geneInfo.csv"

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
source(file.path(lib_path, "GffToGRanges.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# read sample table
sample_table <- 
  readr::read_csv(file = sample_path) %>% 
  dplyr::filter(!(sample_id %in% c("s_ESC_DX_i3_JM7D1.SE", "s_MII_B6_MT_4.SE")))

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read FPKM 
fpkm_df <- 
  readr::read_csv(file.path(outpath, "results", "FPKM.ensembl.GRCm38.89.csv")) %>% 
  dplyr::select(-c(seqnames:gene_description)) %>% 
  dplyr::select_at(.vars = vars(-matches("s_ESC_DX_i3_JM7D1.SE|s_MII_B6_MT_4.SE")))

######################################################## MAIN CODE
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# filter FPKM to only protein coding genes
fpkm_df %<>% 
  dplyr::filter(gene_id %in% protein_genes)

# read summarizedExperiment from RDS file, filter to include only protein coding genes
se <- 
  readRDS(file = file.path(outpath, "Mus_musculus.GRCm38.89.20180305.mESC_oocytes_2018.se.RDS")) %>% 
  .[, -which(str_detect(colnames(.), "s_ESC_DX_i3_JM7D1.SE|s_MII_B6_MT_4.SE"))] %>% 
  .[rownames(.) %in% protein_genes, ]


### diff expression
# set stage to filter and genotype levels
genotype_levels_list <- list(GV = c("ICR_F", "ICR_S"), MII = c("HET", "SOM"))

# get and save results, plot MA in loop
for(filt_stage in c("GV", "MII")){
  
  ### prepare data and variables
  # set genotype levels 
  genotype_levels <- 
    genotype_levels_list[[filt_stage]] %>% 
    sort(.)
  
  # filter sample table
  sample_table_filt <- 
    sample_table %>% 
    dplyr::filter(stage == filt_stage, 
                  genotype %in% genotype_levels) %>% 
    as.data.frame(.) %>% 
    set_rownames(., .$short_name) 
  
  # filter FPKM, get mean across genotype
  fpkm_df_filt <- 
    fpkm_df %>% 
    dplyr::select(which(colnames(.) %in% c("gene_id", sample_table_filt$sample_id))) %>% 
    tidyr::gather(key = sample_id, value = fpkm, -gene_id) %>% 
    dplyr::left_join(., sample_table_filt %>% dplyr::select(sample_id, genotype), by = "sample_id") %>% 
    dplyr::group_by(gene_id, genotype) %>% 
    dplyr::summarise(fpkm = round(mean(fpkm), 3)) %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(key = genotype, value = fpkm) %>% 
    data.table::setnames(., -1, str_c("avgFPKM.", colnames(.)[2:ncol(.)]))   
  
  # filter summarizedExperiment to include only MII, set colData
  se_filt <- se[, colnames(se)[match(basename(sample_table_filt$bam_path), colnames(se))]]
  colData(se_filt) <- DataFrame(sample_table_filt)
  se_filt$genotype <- factor(se_filt$genotype, levels = genotype_levels)
  
  
  ### get differentialy expressed genes
  # run DESeq 
  dds <-
    DESeqDataSet(se_filt, design = ~genotype) %>% 
    DESeq(.)
  
  # shrink results
  dds_shrink <- lfcShrink(dds, coef = 2)
  
  # get results
  results_df <- 
    dds_shrink %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(., var = "gene_id") %>% 
    as.tibble(.) %>% 
    dplyr::arrange(padj) %>% 
    dplyr::left_join(fpkm_df_filt, by = "gene_id") %>% 
    dplyr::left_join(genes_info, by = "gene_id") %>% 
    dplyr::mutate(cooridnates = str_c(seqnames, ":", start, "-", end, "|", strand)) %>% 
    dplyr::select(-c(seqnames:strand)) %>% 
    dplyr::select_at(.vars = vars(gene_id, gene_name, baseMean, log2FoldChange, padj, 
                                  matches("^avgFPKM"), cooridnates, gene_biotype, gene_description)) %T>% 
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.all.csv")))
  
  # write only significant results (padj < 0.1)
  results_df %>% 
    dplyr::filter(padj < 0.1) %>% 
    write_csv(., path = file.path(outpath, "results", str_c("diffExp.", filt_stage, ".", resultsNames(dds)[2], ".GRCm38.89.signif.csv")))
  
  
  ### MA plot
  # create plot
  MA_plot <- 
    results_df %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, padj > 0.1, "no"), 
                  regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
    ggplot(data = ., aes(x = mean, y = lfc, color = regulation, alpha = regulation)) + 
    geom_point(size = 3, shape = 20) +
    scale_x_log10(limits = c(0.1, 1e5), 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    # scale_y_continuous(limits = c(-7, 5),
    #                    breaks = c(-7:5)) +
    scale_colour_manual(values = c(no = "black", down = "#1a75ff", up = "red3")) +
    scale_alpha_manual(values = c(not_sign = 0.3, down = 1, up = 1)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    xlab("Mean expression") +
    ylab(str_c("log2FoldChange ", resultsNames(dds)[2])) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
          axis.title.y = element_text(size = 15, vjust = 0.3), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, "results", str_c("MAplot.", filt_stage, ".", resultsNames(dds)[2], ".png")),
         plot = MA_plot, width = 15, height = 10)
  
  
  ### interactive MA plot
  # data for plot
  interactive_plot_df <-
    results_df %>%
    dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, 
                                  matches("^avgFPKM"))) %>%
    dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                  regulation = ifelse(lfc > 0, "up", "down"),
                  regulation = replace(regulation, padj > 0.1, "no"),
                  regulation = factor(regulation, levels = c("no", "down", "up")),
                  gene_description = str_replace(gene_description, " \\[.*", ""), 
                  gene_description = replace(gene_description, is.na(gene_description), "")) %>%
    dplyr::arrange(regulation)
  
  # make interactive MA plot
  p <-
    plotly::plot_ly(data = interactive_plot_df,
                    x = ~mean,
                    y = ~lfc,
                    text = ~paste("</br> log2FC: ", round(lfc, 3), 
                                  "</br> mean exp.: ", round(mean, 3), 
                                  str_c("</br> ", genotype_levels[1], " avg. FPKM:"), get(str_c("avgFPKM.", genotype_levels[1])),
                                  str_c("</br> ", genotype_levels[2], " avg. FPKM:"), get(str_c("avgFPKM.", genotype_levels[2])),
                                  "</br>", gene_id,
                                  "</br>", gene_name,
                                  "</br>", gene_description),
                    color = ~regulation,
                    colors = c("gray32", "#1a75ff", "red3"),
                    alpha = 0.75,
                    hoverinfo = "text") %>%
    add_markers() %>%
    layout(xaxis = list(title = "mean expression", type = "log"),
           yaxis = list(title = str_c("log2FoldChange ", resultsNames(dds)[2])))
  
  # save as html widget
  htmlwidgets::saveWidget(as_widget(p),
                          file = file.path(outpath, "results", str_c("MAplot.", filt_stage, ".", resultsNames(dds)[2], ".html")),
                          selfcontained = T)
  
}

