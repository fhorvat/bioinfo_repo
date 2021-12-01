### INFO: Do the differential expression analysis for mESC and oocytes sequenced in February 2018
### DATE: Fri Nov 30 18:24:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Analysis/mESC_DX")

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
library(openxlsx)
library(VennDiagram)
library(geneplotter)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set input path
inpath <- getwd()

# set output path
outpath <- getwd()

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/mESC_oocytes_2018/Data/Documentation/mESC_oocytes_2018.sample_table.csv"

# genes info
genes_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.geneInfo.csv"

# summarizedExperiment path
se_path <- file.path(inpath, "mESC_DX.GRCm38.91.reducedExons.summarizedOverlaps.RDS")

# FPKM path
fpkm_path <- file.path(inpath, "mESC_DX.GRCm38.91.reducedExons.FPKM.csv")

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path) 

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read summarizedExperiment
se <- readRDS(file = se_path) 

# read FPKM 
fpkm_df <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
#### PREPARE DATA ####
# mutate gene info
genes_info %<>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ")

# filter sample table
sample_table_dds <- 
  sample_table %>% 
  dplyr::filter(str_detect(sample_id, "^s_ESC_DX_.*")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE"), 
                genotype_simple = str_extract(sample_id, "JME|JM7D|JM71|R1BD"), 
                genotype_very_simple = ifelse(str_detect(genotype, "JME"), "WT", "DcrX")) %>% 
  dplyr::mutate(genotype_simple = factor(genotype_simple, levels = c("JME", "JM7D", "JM71", "R1BD")), 
                genotype_very_simple = factor(genotype_very_simple, levels = c("WT", "DcrX"))) %>% 
  dplyr::filter(sample_id != "s_ESC_DX_i3_JM7D1") %>%
  as.data.frame(.) %>% 
  set_rownames(., .$sample_id) 

# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# get long FPKM table
fpkm_long <- 
  fpkm_df %>% 
  dplyr::filter(gene_id %in% protein_genes) %>% 
  dplyr::select(-c(s_ESC_DX_i3_JM7D1, coordinates:gene_description)) %>% 
  tidyr::gather(key = sample_id, value = fpkm, -gene_id) %>% 
  dplyr::left_join(., sample_table_dds %>% dplyr::select(sample_id, genotype, genotype_simple, genotype_very_simple), by = "sample_id") %>% 
  dplyr::mutate(fpkm = replace(fpkm, is.na(fpkm), 0)) 

# get mean FPKM - genotype simple
fpkm_mean.genotype_simple <- 
  fpkm_long %>% 
  dplyr::group_by(gene_id, genotype_simple) %>% 
  dplyr::summarise(fpkm = round(mean(fpkm), 3)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = genotype_simple, value = fpkm) %>% 
  data.table::setnames(., -1, str_c(colnames(.)[2:ncol(.)], ".avgFPKM"))   

# get mean FPKM - genotype very simple
fpkm_mean.genotype_very_simple <- 
  fpkm_long %>% 
  dplyr::group_by(gene_id, genotype_very_simple) %>% 
  dplyr::summarise(fpkm = round(mean(fpkm), 3)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = genotype_very_simple, value = fpkm) %>% 
  data.table::setnames(., -1, str_c(colnames(.)[2:ncol(.)], ".avgFPKM"))   

### summarizedExperiment
# change colnames
colnames(se) <- str_remove(colnames(se), "\\.SE\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam")

# read summarizedExperiment from RDS file, filter to include only protein coding genes
se_filt <- 
  se %>% 
  .[, !(str_detect(colnames(.), "s_ESC_DX_i3_JM7D1"))] %>%
  .[rownames(.) %in% protein_genes, ] %>% 
  .[, colnames(.)[match(rownames(sample_table_dds), colnames(.))]]

# add column data to SE
colData(se_filt) <- DataFrame(sample_table_dds)


#### DESEQ2 DIFFERENTIAL EXPRESSION - GENOTYPE SIMPLE ####
# set 3 simple genotypes as vector
genotypes_simple <- c("JM7D", "JM71", "R1BD")

# make DESeqDataSet
dds <- 
  DESeqDataSet(se_filt, design = ~genotype_simple) %>% 
  DESeq(.)

# get results for 3 different DcrX genotypes compared to WT
results_all <- purrr::map(genotypes_simple, function(genotype){
  
  # shrink results
  dds_shrink <- lfcShrink(dds, contrast = c("genotype_simple", genotype, "JME"))
  
  # get results
  results_df <- 
    dds_shrink %>% 
    as.data.frame(.) %>% 
    as.tibble(., rownames = "gene_id") %>% 
    dplyr::arrange(padj) %>% 
    dplyr::left_join(fpkm_mean.genotype_simple %>% dplyr::select_at(vars("gene_id", contains("JME"), contains(genotype))), by = "gene_id") %>% 
    dplyr::left_join(genes_info, by = "gene_id") %>% 
    dplyr::select_at(vars(gene_id, gene_name, baseMean, log2FoldChange, padj, 
                          matches(".*avgFPKM"), coordinates, gene_biotype, gene_description))
  
}) %>% 
  magrittr::set_names(genotypes_simple)

# get significant results
results_significant <- purrr::map(results_all, function(results_df){
  
  # filter by p-adjusted (padj < 0.1)
  results_df_signif <- 
    results_df %>%
    dplyr::filter(padj <= 0.1)
  
})


### save results
# save all results
write.xlsx(x = results_all, 
           file = file.path(outpath, str_c("diffExp.genotype_simple.GRCm38.91.all.xlsx")), 
           asTable = rep(TRUE, length(results_all)),
           sheetName = str_c(genotypes_simple, "_vs_JME"), 
           colWidths = "auto")

# save significant results
write.xlsx(x = results_significant, 
           file = file.path(outpath, str_c("diffExp.genotype_simple.GRCm38.91.significant.xlsx")), 
           asTable = rep(TRUE, length(results_significant)),
           sheetName = str_c(genotypes_simple, "_vs_JME"), 
           colWidths = "auto")


### MA plots
# save as .png
purrr::map(names(results_all), function(genotype){
  
  # create plot table
  MA_df <- 
    results_all[[genotype]] %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, padj > 0.1, "no"), 
                  regulation = factor(regulation, levels = c("no", "down", "up")))
  
  # plot 
  MA_plot <- 
    ggplot(data = MA_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation)) + 
    geom_point(size = 3, shape = 20) +
    scale_x_log10(limits = c(0.1, 1e5), 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_continuous(limits = c(-4, 4), breaks = c(-4:4)) +
    scale_colour_manual(values = c(no = "gray32", down = "#1a75ff", up = "red3")) +
    scale_alpha_manual(values = c(no = 1, down = 1, up = 1)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    xlab("Mean expression") +
    ylab(str_c("log2FoldChange ", genotype, " vs JME")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
          axis.title.y = element_text(size = 15, vjust = 0.3), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("MAplot.genotype_simple.", genotype, "_vs_JME.png")),
         plot = MA_plot, width = 15, height = 10)
  
})

# save as .html
purrr::map(names(results_all), function(genotype){
  
  # create plot table
  MA_df <- 
    results_all[[genotype]] %>% 
    dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, matches("\\.avgFPKM$"))) %>%
    dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                  regulation = ifelse(lfc > 0, "up", "down"),
                  regulation = replace(regulation, padj > 0.1, "no"),
                  regulation = factor(regulation, levels = c("no", "down", "up")),
                  gene_description = str_replace(gene_description, " \\[.*", ""),
                  gene_description = replace(gene_description, is.na(gene_description), "")) %>%
    dplyr::arrange(regulation)
  
  # make interactive MA plot
  if(nrow(MA_df) > 0){
    
    interactive_ma_plot <-
      plotly::plot_ly(data = MA_df,
                      x = ~mean,
                      y = ~lfc,
                      text = ~paste("</br> log2FC: ", round(lfc, 3),
                                    "</br> mean exp.: ", round(mean, 3),
                                    str_c("</br> ", genotype, " avg. FPKM:"), get(str_c(genotype, ".avgFPKM")),
                                    str_c("</br> ", "JME", " avg. FPKM:"), get(str_c("JME", ".avgFPKM")),
                                    "</br>", gene_id,
                                    "</br>", gene_name,
                                    "</br>", gene_description),
                      color = ~regulation,
                      colors = c("gray32", "#1a75ff", "red3"),
                      alpha = 0.75,
                      hoverinfo = "text") %>%
      add_markers() %>%
      layout(xaxis = list(title = "mean expression", type = "log"),
             yaxis = list(title = str_c("log2FoldChange ", genotype, " vs. JME")))
    
    # save as html widget
    htmlwidgets::saveWidget(as_widget(interactive_ma_plot),
                            file = file.path(outpath, str_c("MAplot.genotype_simple.", genotype, "_vs_JME.html")),
                            selfcontained = T)
    
  }
  
})


### scatter plots
# save as .png
purrr::map(names(results_all), function(genotype){
  
  # create plot table
  scatter_df <- 
    results_all[[genotype]] %>% 
    dplyr::select(WT = JME.avgFPKM, DcrX = str_c(genotype, ".avgFPKM"), lfc = log2FoldChange, padj, gene_id) %>% 
    dplyr::mutate_at(vars("WT", "DcrX"), funs(log2(.))) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, padj > 0.1, "no"), 
                  regulation = factor(regulation, levels = c("no", "down", "up"))) %>% 
    dplyr::arrange(regulation)
  
  # plot 
  scatter_plot <- 
    ggplot(data = scatter_df, aes(x = WT, y = DcrX, color = regulation, alpha = regulation)) + 
    geom_point(size = 3, shape = 20) +
    # scale_x_log10(limits = c(0.1, 1e5), 
    #               breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    # scale_y_continuous(limits = c(-4, 4), breaks = c(-4:4)) +
    scale_colour_manual(values = c(no = "gray32", down = "#1a75ff", up = "red3")) +
    scale_alpha_manual(values = c(no = 1, down = 1, up = 1)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    xlab("log2 FPKM in JME") +
    ylab(str_c("log2 FPKM in ", genotype)) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
          axis.title.y = element_text(size = 15, vjust = 0.3), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("scatterPlot.genotype_simple.", genotype, "_vs_JME.png")),
         plot = scatter_plot, width = 10, height = 10)
  
})


#### Venn diagram ####
# down-regulated genes
down_genes <- purrr::map(results_significant, function(x){
  x %>% dplyr::filter(log2FoldChange < 0) %$% gene_id
})

# up-regulated genes
up_genes <- lapply(results_significant, function(x){
  x %>% dplyr::filter(log2FoldChange > 0) %$% gene_id
})

# down-regulated 
png(file = file.path(outpath, "Venn.genotype_simple.significant_downregulated.DcrX_vs_WT.GRCm38.91.png"), 
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
png(file = file.path(outpath, "Venn.genotype_simple.significant_upregulated.DcrX_vs_WT.GRCm38.91.png"), 
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








# #### DESEQ2 DIFFERENTIAL EXPRESSION - GENOTYPE VERY SIMPLE ####
# # set genotype
# genotype_very_simple <- "DcrX"
# 
# # make DESeqDataSet
# dds <- 
#   DESeqDataSet(se_filt, design = ~genotype_very_simple) %>% 
#   DESeq(.)
# 
# # get results for DcrX genotypes compared to WT
# results_all <- purrr::map(genotype_very_simple, function(genotype){
#   
#   # shrink results
#   dds_shrink <- lfcShrink(dds, contrast = c("genotype_very_simple", genotype_very_simple, "WT"))
#   
#   # get results
#   results_df <- 
#     dds_shrink %>% 
#     as.data.frame(.) %>% 
#     as.tibble(., rownames = "gene_id") %>% 
#     dplyr::arrange(padj) %>% 
#     dplyr::left_join(fpkm_mean.genotype_very_simple %>% dplyr::select_at(vars("gene_id", contains("WT"), contains(genotype_very_simple))), by = "gene_id") %>% 
#     dplyr::left_join(genes_info, by = "gene_id") %>% 
#     dplyr::select_at(vars(gene_id, gene_name, baseMean, log2FoldChange, padj, 
#                           matches(".*avgFPKM"), coordinates, gene_biotype, gene_description))
#   
# }) %>% 
#   magrittr::set_names(genotype_very_simple)
# 
# # get significant results
# results_significant <- purrr::map(results_all, function(results_df){
#   
#   # filter by p-adjusted (padj < 0.1)
#   results_df_signif <- 
#     results_df %>%
#     dplyr::filter(padj <= 0.1)
#   
# })
# 
# 
# ### save results
# # save all results
# write.xlsx(x = results_all, 
#            file = file.path(outpath, str_c("diffExp.genotype_very_simple.GRCm38.91.all.xlsx")), 
#            asTable = rep(TRUE, length(results_all)),
#            sheetName = str_c(genotype_very_simple, "_vs_WT"), 
#            colWidths = "auto")
# 
# # save significant results
# write.xlsx(x = results_significant, 
#            file = file.path(outpath, str_c("diffExp.genotype_very_simple.GRCm38.91.significant.xlsx")), 
#            asTable = rep(TRUE, length(results_significant)),
#            sheetName = str_c(genotype_very_simple, "_vs_WT"), 
#            colWidths = "auto")
# 
# 
# ### MA plots
# # save as .png
# purrr::map(names(results_all), function(genotype){
#   
#   # create plot table
#   MA_df <- 
#     results_all[[genotype]] %>% 
#     dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
#     dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
#                   regulation = ifelse(lfc > 0, "up", "down"), 
#                   regulation = replace(regulation, padj > 0.1, "no"), 
#                   regulation = factor(regulation, levels = c("no", "down", "up")))
#   
#   # plot 
#   MA_plot <- 
#     ggplot(data = MA_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation)) + 
#     geom_point(size = 3, shape = 20) +
#     scale_x_log10(limits = c(0.1, 1e5), 
#                   breaks = scales::trans_breaks("log10", function(x) 10^x),
#                   labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
#     # scale_y_continuous(limits = c(-7, 5), breaks = c(-7:5)) +
#     scale_colour_manual(values = c(no = "gray32", down = "#1a75ff", up = "red3")) +
#     scale_alpha_manual(values = c(no = 1, down = 1, up = 1)) +
#     guides(color = guide_legend(override.aes = list(size = 5))) +
#     xlab("Mean expression") +
#     ylab(str_c("log2FoldChange ", resultsNames(dds)[2])) +
#     theme_bw() +
#     theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
#           axis.title.y = element_text(size = 15, vjust = 0.3), 
#           axis.text.x = element_text(size = 15), 
#           axis.text.y = element_text(size = 15),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())
#   
#   # save plot
#   ggsave(filename = file.path(outpath, str_c("MAplot.genotype_very_simple.", genotype_very_simple, "_vs_WT.png")),
#          plot = MA_plot, width = 15, height = 10)
#   
# })
# 
# 
# ### interactive MA plots
# # save as .html
# purrr::map(names(results_all), function(genotype){
#   
#   # create plot table
#   MA_df <- 
#     results_all[[genotype]] %>% 
#     dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, matches("\\.avgFPKM$"))) %>%
#     dplyr::mutate(padj = replace(padj, is.na(padj), 1),
#                   regulation = ifelse(lfc > 0, "up", "down"),
#                   regulation = replace(regulation, padj > 0.1, "no"),
#                   regulation = factor(regulation, levels = c("no", "down", "up")),
#                   gene_description = str_replace(gene_description, " \\[.*", ""),
#                   gene_description = replace(gene_description, is.na(gene_description), "")) %>%
#     dplyr::arrange(regulation)
#   
#   # make interactive MA plot
#   if(nrow(MA_df) > 0){
#     
#     interactive_ma_plot <-
#       plotly::plot_ly(data = MA_df,
#                       x = ~mean,
#                       y = ~lfc,
#                       text = ~paste("</br> log2FC: ", round(lfc, 3),
#                                     "</br> mean exp.: ", round(mean, 3),
#                                     str_c("</br> ", genotype, " avg. FPKM:"), get(str_c(genotype, ".avgFPKM")),
#                                     str_c("</br> ", "WT", " avg. FPKM:"), get(str_c("WT", ".avgFPKM")),
#                                     "</br>", gene_id,
#                                     "</br>", gene_name,
#                                     "</br>", gene_description),
#                       color = ~regulation,
#                       colors = c("gray32", "#1a75ff", "red3"),
#                       alpha = 0.75,
#                       hoverinfo = "text") %>%
#       add_markers() %>%
#       layout(xaxis = list(title = "mean expression", type = "log"),
#              yaxis = list(title = str_c("log2FoldChange ", genotype, " vs. WT")))
#     
#     # save as html widget
#     htmlwidgets::saveWidget(as_widget(interactive_ma_plot),
#                             file = file.path(outpath, str_c("MAplot.genotype_very_simple.", genotype_very_simple, "_vs_WT.html")),
#                             selfcontained = T)
#     
#   }
#   
# })
