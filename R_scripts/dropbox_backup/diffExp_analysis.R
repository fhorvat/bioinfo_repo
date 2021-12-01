#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: get FPKM and count tables
### DATE: Wed Jun 19 16:48:28 2019
### AUTHOR: Filip Horvat
options(bitmapType = "cairo")

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
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# list of read stats and mapped logs
se_path <- args$se_path
fpkm_path <- args$fpkm_path
sample_table_path <- args$sample_table_path
genes_info_path <- args$genes_info_path
grouping_variables <- args$grouping_variables
results_groups <- args$results_groups

######################################################## READ DATA
# summarizedExperiment 
se <- readRDS(se_path)

# read mean FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

# read sample table
sample_table <- readr::read_csv(sample_table_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
### prepare data
# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# create vector of plotly symbols in ggplot shape order
ploty_symbols <- c("square", "circle", "cross", "x", "diamond", "square-open", "circle-open", "diamond-open")


### exploratory analysis 
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  as.data.table(.) %>% 
  .[, c("sample_id", grouping_variables), with = F] %>% 
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage and protein coding genes
se_filt <- se[rownames(se) %in% protein_genes, ]
colnames(se_filt) <- str_remove(colnames(se_filt), "\\.bam$")
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
# create stargin plot object
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
ggsave(filename = file.path(outpath, str_c("diffExp", str_c(grouping_variables, collapse = "_"), "rlog.PC1_2.PCA_plot.png", sep = ".")),
       plot = pca_plot, width = 10, height = 10)

# add labels
pca_plot <- 
  pca_plot +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save labeled plot
ggsave(filename = file.path(outpath, str_c("diffExp", str_c(grouping_variables, collapse = "_"), "rlog.PC1_2.labeled.PCA_plot.png", sep = ".")),
       plot = pca_plot, width = 10, height = 10)


### 3D PCA plot
# data for plot
pca_3Dplot <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         PC3 = pca$x[, 3],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))

# if there is only one grouping variable use only color, if there is more use also a shape
if(length(grouping_variables) == 1){
  
  # make interactive 3D PCA plot
  p <-
    plotly::plot_ly(data = pca_3Dplot,
                    x = ~PC1,
                    y = ~PC2,
                    z = ~PC3,
                    color = ~pull(pca_3Dplot[, grouping_variables[1]]),
                    colors = gg_color_hue(pull(pca_3Dplot[, grouping_variables[1]])%>% unique(.) %>% length(.)), 
                    mode = "markers", 
                    marker = list(size = 10)
                    
    ) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = str_c("PC1: ", round(percentVar[1] * 100), "% variance"), range = c(-40, 40)),
                        yaxis = list(title = str_c("PC2: ", round(percentVar[2] * 100), "% variance"), range = c(-40, 40)), 
                        zaxis = list(title = str_c("PC3: ", round(percentVar[3] * 100), "% variance"), range = c(-40, 40))))
  
  # save as html widget
  htmlwidgets::saveWidget(plotly::as_widget(p),
                          file = file.path(outpath, str_c("diffExp", str_c(grouping_variables, collapse = "_"), "rlog.PC1_2.PCA_plot.html", sep = ".")),
                          selfcontained = T)
}else{
  
  # make interactive 3D PCA plot
  p <-
    plotly::plot_ly(data = pca_3Dplot,
                    x = ~PC1,
                    y = ~PC2,
                    z = ~PC3,
                    color = ~pull(pca_3Dplot[, grouping_variables[1]]),
                    colors = gg_color_hue(pull(pca_3Dplot[, grouping_variables[1]])%>% unique(.) %>% length(.)), 
                    symbol = ~pull(pca_3Dplot[, grouping_variables[2]]), 
                    symbols = ploty_symbols[1:(pull(pca_3Dplot[, grouping_variables[2]]) %>% unique(.) %>% length(.))],
                    mode = "markers", 
                    marker = list(size = 10)
                    
    ) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = str_c("PC1: ", round(percentVar[1] * 100), "% variance"), range = c(-40, 40)),
                        yaxis = list(title = str_c("PC2: ", round(percentVar[2] * 100), "% variance"), range = c(-40, 40)), 
                        zaxis = list(title = str_c("PC3: ", round(percentVar[3] * 100), "% variance"), range = c(-40, 40))))
  
  # save as html widget
  htmlwidgets::saveWidget(plotly::as_widget(p),
                          file = file.path(outpath, str_c("diffExp", str_c(grouping_variables, collapse = "_"), "rlog.PC1_2.PCA_plot.html", sep = ".")),
                          selfcontained = T)
}


### distance heatmap
# calculate distance
dist_df <-
  rlog_df %>%
  t(.) %>%
  dist(.)

# make matrix
dist_matrix <- as.matrix(dist_df)
colnames(dist_matrix) <- NULL

# plot
pheatmap::pheatmap(dist_matrix,
                   clustering_distance_rows = dist_df,
                   clustering_distance_cols = dist_df,
                   col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   filename = file.path(outpath,  str_c("diffExp", str_c(grouping_variables, collapse = "_"), "rlog.distance_heatmap.png", sep = ".")),
                   height = 10,
                   width = 12)


### DESeq2
# run DESeq
dds_deseq <- DESeq(dds)

### loop through results from results_groups variable
# check whether to do diff. exp. analysis
if(!is.null(results_groups)){
  
  # create results workbooks
  wb_all <- openxlsx::createWorkbook()
  wb_significant <- openxlsx::createWorkbook()
  
  # loop
  for(result in results_groups){
    
    # shape result
    result_clean <- 
      str_split(result, pattern = ",") %>% 
      unlist(.)
    
    # check
    if((length(result_clean) == 2) & (all(result_clean %in% sample_table_dds$grouped_variables))){
      
      # get results, shrink logFC
      dds_shrink <- lfcShrink(dds_deseq, contrast = c("grouped_variables", result_clean[1], result_clean[2]))
      
      # add sheet to workbook
      openxlsx::addWorksheet(wb = wb_all, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]))
      
      # get results table, add to sheet
      results_df <-
        dds_shrink %>%
        as_tibble(., rownames = "gene_id") %>%
        dplyr::arrange(padj) %>%
        dplyr::left_join(fpkm_tb %>% dplyr::select(-one_of(sample_table_dds %>% dplyr::filter(!(grouped_variables %in% result_clean)) %$% grouped_variables %>% unique(.))), 
                         by = "gene_id") %>%
        dplyr::mutate(comparison = str_c(result_clean[1], "_vs_", result_clean[2])) %>% 
        setnames(., old = result_clean, new = str_c(result_clean, ".FPKM")) %T>% 
        openxlsx::writeData(wb = wb_all, 
                            sheet = str_c(result_clean[1], "_vs_", result_clean[2]),
                            x = .)
      
      ### write only significant results, padj < 0.1
      # filter table
      results_df_sign <- 
        results_df %>% 
        dplyr::filter(padj < 0.1)
      
      # check and write
      if(nrow(results_df_sign) > 0){
        
        # open sheet
        openxlsx::addWorksheet(wb = wb_significant, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]))
        
        # write data
        openxlsx::writeData(wb = wb_significant, 
                            sheet = str_c(result_clean[1], "_vs_", result_clean[2]),
                            x = results_df_sign)
        
      }
      
      
      ### MA plot - ggplot2
      # data for plot
      plot_df <-
        results_df %>%
        dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      sign = ifelse(padj < 0.1, "yes", "no"),
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
        scale_y_continuous(limits = c(-7, 7),
                           breaks = c(-7:7)) +
        scale_colour_manual(values = c(not_sign = "gray30", up = "red3", down = "blue3")) +
        guides(color = FALSE) +
        xlab("mean expression") +
        ylab("log2FC") +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      # save plot
      ggsave(filename = file.path(outpath, str_c("diffExp", 
                                                 str_c(grouping_variables, collapse = "_"), 
                                                 str_c(result_clean[1], "_vs_", result_clean[2]), 
                                                 "DESeq2.MA_plot.png", sep = ".")),
             plot = ma_plot, width = 10, height = 10)
      
      
      ### MA plot - interactive plot.ly
      # data for plot
      plot_df <-
        results_df %>%
        dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, contains("FPKM"))) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      regulation = ifelse(lfc > 0, "up", "down"),
                      regulation = replace(regulation, padj > 0.1, "no"),
                      regulation = factor(regulation, levels = c("no", "down", "up")),
                      gene_description = str_remove(gene_description, " \\[.*"),
                      gene_description = replace(gene_description, is.na(gene_description), "")) %>%
        dplyr::arrange(regulation)
      
      # make interactive MA plot
      if(nrow(plot_df) > 0){
        
        interactive_ma_plot <-
          plotly::plot_ly(data = plot_df,
                          x = ~mean,
                          y = ~lfc,
                          text = ~paste("</br> log2FC: ", round(lfc, 3),
                                        "</br> mean exp.: ", round(mean, 3),
                                        str_c("</br> ", result_clean[1], " FPKM:"), get(str_c(result_clean[1], ".FPKM")),
                                        str_c("</br> ", result_clean[2], " FPKM:"), get(str_c(result_clean[2], ".FPKM")),
                                        "</br>", gene_id,
                                        "</br>", gene_name,
                                        "</br>", gene_description),
                          color = ~regulation,
                          colors = c("gray32", "#1a75ff", "red3"),
                          alpha = 0.75,
                          hoverinfo = "text") %>%
          add_markers() %>%
          layout(xaxis = list(title = "mean expression", type = "log"),
                 yaxis = list(title = "log2FoldChange"))
        
        # save as html widget
        htmlwidgets::saveWidget(plotly::as_widget(interactive_ma_plot),
                                file = file.path(outpath, str_c("diffExp", 
                                                                str_c(grouping_variables, collapse = "_"), 
                                                                str_c(result_clean[1], "_vs_", result_clean[2]), 
                                                                "DESeq2.MA_plot.html", sep = ".")), 
                                selfcontained = T)
        
      }
      
      
      ### Vulcano plot - ggplot2
      # set p-adjusted and log2FC cutoff
      padj_cut <- 0.1
      lfc_cut <- 1.5
      
      # data for plot
      plot_df <-
        results_df %>%
        dplyr::select(lfc = log2FoldChange, padj, gene_id, gene_name) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      padj_sign = ifelse(padj <= padj_cut, T, F),
                      lfc_sign = ifelse((abs(lfc) > lfc_cut), T, F)) %>%
        dplyr::mutate(regulation = "not_sign",
                      regulation = replace(regulation, lfc_sign, "fold_change"),
                      regulation = replace(regulation, padj_sign, "p_value"),
                      regulation = replace(regulation, lfc_sign & padj_sign, "p_value_fold_change"),
                      regulation = factor(regulation, levels = c("not_sign", "fold_change", "p_value", "p_value_fold_change"))) %>%
        dplyr::arrange(regulation)
      
      # plot
      vulcano_plot <-
        ggplot(plot_df, aes(x = lfc, y = -log10(padj), color = regulation)) +
        geom_point(size = 1.5, shape = 19, alpha = 0.5) +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "longdash", colour = "black", size = 0.4) +
        geom_hline(yintercept = -log10(padj_cut), linetype = "longdash", colour = "black", size = 0.4) +
        geom_text(data = plot_df %>% dplyr::filter(padj_sign & lfc_sign),
                  aes(label = plot_df %>% dplyr::filter(padj_sign & lfc_sign) %$% gene_name),
                  check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
                  colour = "black", fontface = "plain") +
        scale_colour_manual(values = c(not_sign = "gray30", fold_change = "forestgreen", p_value = "royalblue", p_value_fold_change = "red2")) +
        scale_x_continuous(limits = c(-8, 8)) +
        guides(color = F, alpha = F) +
        xlab("log2FC") +
        ylab("-log10 p-adjusted value") +
        theme_bw(base_size = 24) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        theme(legend.position = "top")
      
      # save plot
      ggsave(filename = file.path(outpath, str_c("diffExp", 
                                                 str_c(grouping_variables, collapse = "_"), 
                                                 str_c(result_clean[1], "_vs_", result_clean[2]), 
                                                 "DESeq2.vulcano_plot.png", sep = ".")),
             plot = vulcano_plot, width = 10, height = 10)
      
      
      ### Vulcano plot - interactive plot.ly
      # data for plot
      plot_df <-
        results_df %>%
        dplyr::select_at(.vars = vars(lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, contains("FPKM"))) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      padj_sign = ifelse(padj <= padj_cut, T, F),
                      lfc_sign = ifelse((abs(lfc) > lfc_cut), T, F)) %>%
        dplyr::mutate(regulation = "not_sign",
                      regulation = replace(regulation, lfc_sign, "fold_change"),
                      regulation = replace(regulation, padj_sign, "p_value"),
                      regulation = replace(regulation, lfc_sign & padj_sign, "p_value_fold_change"),
                      regulation = factor(regulation, levels = c("not_sign", "fold_change", "p_value", "p_value_fold_change"))) %>%
        dplyr::mutate(gene_description = str_remove(gene_description, " \\[.*"),
                      gene_description = replace(gene_description, is.na(gene_description), "")) %>%
        dplyr::arrange(regulation)
      
      
      # make interactive Vulcano plot
      if(nrow(plot_df) > 0){
        
        interactive_vulcano_plot <-
          plotly::plot_ly(data = plot_df,
                          x = ~lfc,
                          y = ~-log10(padj),
                          text = ~paste("</br> log2FC: ", round(lfc, 3),
                                        "</br> p-adjusted: ", round(padj, 3),
                                        str_c("</br> ", result_clean[1], " FPKM:"), get(str_c(result_clean[1], ".FPKM")),
                                        str_c("</br> ", result_clean[2], " FPKM:"), get(str_c(result_clean[2], ".FPKM")),
                                        "</br>", gene_id,
                                        "</br>", gene_name,
                                        "</br>", gene_description),
                          color = ~regulation,
                          colors = c("gray30", "forestgreen", "royalblue", "red2"),
                          alpha = 0.75,
                          hoverinfo = "text") %>%
          add_markers() %>%
          layout(xaxis = list(title = "log2FC"),
                 yaxis = list(title = str_c("-log10 p-adjusted value")))
        
        # save as html widget
        htmlwidgets::saveWidget(plotly::as_widget(interactive_vulcano_plot),
                                file = file.path(outpath, str_c("diffExp", 
                                                                str_c(grouping_variables, collapse = "_"), 
                                                                str_c(result_clean[1], "_vs_", result_clean[2]), 
                                                                "DESeq2.vulcano_plot.html", sep = ".")),
                                selfcontained = T)
        
      }
      
    }
    
  }
  
  # save workbooks to outdir
  openxlsx::saveWorkbook(wb = wb_all, 
                         file = file.path(outpath, str_c("diffExp", str_c(grouping_variables, collapse = "_"), "DESeq2.all_results.xlsx", sep = ".")), 
                         overwrite = TRUE)
  
  openxlsx::saveWorkbook(wb = wb_significant, 
                         file = file.path(outpath, str_c("diffExp", str_c(grouping_variables, collapse = "_"), "DESeq2.significant_results.xlsx", sep = ".")), 
                         overwrite = TRUE)
  
}


