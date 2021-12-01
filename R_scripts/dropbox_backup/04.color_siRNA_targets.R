### INFO: 
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Ago2_KO/datasets/2019_Oct/Analysis/expression/siRNA_targets")

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
### in and out
# set inpath
inpath <- file.path(getwd(), "../results.ensembl.93.GRCm38.p6.20180919.UCSCseqnames")

# set outpath
outpath <- getwd()

# siRNA targets list
sirna_targets_path <- file.path(outpath, "2019 RNAi &  piRNA file 2 siRNA target list.csv")

# results path
results_path <- file.path(inpath, "all_biotype.diffExp.DESeq2.genotype.all_results.xlsx")

######################################################## READ DATA
# read siRNA targets table
sirna_targets <- readr::read_csv(sirna_targets_path)

# read differential expression results
results_list <- purrr::map(openxlsx::getSheetNames(results_path), function(sheet){
  
  # read sheet by sheet, convert to tibble
  openxlsx::read.xlsx(results_path, sheet) %>% 
    as_tibble(.)
  
}) %>% 
  magrittr::set_names(., openxlsx::getSheetNames(results_path))

######################################################## MAIN CODE
### create static MA plots
# get axis limits
results_limits <- 
  results_list %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::summarise(x_limit = baseMean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.), 
                   y_limit = log2FoldChange %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# loop through results, plot and save
invisible(purrr::map(names(results_list), function(result){
  
  ## prepare results
  # get results table
  results_df <- results_list[[result]]
  
  # shape result
  result_clean <- 
    str_split(result, pattern = "_vs_") %>% 
    unlist(.)
  
  ## MA plot
  # data for plot
  plot_df <-
    results_df %>%
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>%
    dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                  padj = replace(padj, padj == 0, .Machine$double.xmin),
                  sign = ifelse(padj < 0.1, "yes", "no"),
                  regulation = ifelse(lfc > 0, "up", "down"),
                  regulation = replace(regulation, padj > 0.1, "no"),
                  regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
    dplyr::mutate(sirna_target = ifelse(gene_id %in% sirna_targets$gene_id, "yes", "no"), 
                  sirna_target = factor(sirna_target, levels = c("no", "yes"))) %>% 
    dplyr::arrange(sirna_target)
  
  # plot
  ma_plot <-
    ggplot() + 
    geom_point(data = plot_df, aes(x = mean, y = lfc, color = sirna_target, alpha = sirna_target), size = 2.5, shape = 20) +
    scale_x_log10(limits = c(0.01, results_limits$x_limit),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit), 
                       breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
    scale_colour_manual(labels = c(no = "not target", yes = "siRNA target"), 
                        values = c(no = "gray50", yes = "red2")) +
    scale_alpha_manual(values = c(no = 0.5, yes = 1)) +      
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("gray50", "red2"))), 
           alpha = F) +
    xlab("mean expression") +
    ylab(str_c("log2 fold change: ", result_clean[1], " / ", result_clean[2], "\n") %>% str_replace_all(., "_", " ")) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.title.x = element_text(size = 13), 
          axis.title.y = element_text(size = 13)) +
    theme(legend.title = element_blank())
  
  # turns off axis titles and legend
  ma_plot <-
    ma_plot +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, 
                              str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), 
                                    "plot", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                    str_c(result_clean[1], "_vs_", result_clean[2]), 
                                    "siRNA_targets", 
                                    "png", sep = ".")),
         plot = ma_plot, width = 12, height = 10)
  
}))



### create interactive MA plots
# loop through results
invisible(purrr::map(names(results_list), function(result){
  
  ## prepare results
  # get results table
  results_df <- results_list[[result]]
  
  # shape result
  result_clean <- 
    str_split(result, pattern = "_vs_") %>% 
    unlist(.)
  
  
  ### MA plot - interactive plot.ly
  # data for plot
  plot_df <-
    results_df %>%
    dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_description, contains("FPKM"))) %>%
    dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                  padj = replace(padj, padj == 0, .Machine$double.xmin),
                  regulation = ifelse(lfc > 0, "up", "down"),
                  regulation = replace(regulation, padj > 0.1, "no"),
                  regulation = factor(regulation, levels = c("no", "down", "up")),
                  gene_description = str_remove(gene_description, " \\[.*"),
                  gene_description = replace(gene_description, is.na(gene_description), "")) %>%
    dplyr::mutate(sirna_target = ifelse(gene_id %in% sirna_targets$gene_id, "yes", "no"), 
                  sirna_target = factor(sirna_target, levels = c("no", "yes"))) %>% 
    dplyr::arrange(sirna_target)
  
  # plot
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
                                    "</br>", gene_description),
                      color = ~regulation,
                      colors = c("gray32", "red3"),
                      alpha = 0.75,
                      hoverinfo = "text") %>%
      add_markers() %>%
      layout(xaxis = list(title = "mean expression", type = "log"),
             yaxis = list(title = "log2FoldChange"))
    
    # save as html widget
    htmlwidgets::saveWidget(plotly::as_widget(interactive_ma_plot),
                            file = file.path(outpath, 
                                             str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), 
                                                   "interactive", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                   str_c(result_clean[1], "_vs_", result_clean[2]),
                                                   "siRNA_targets", 
                                                   "html", sep = ".")), 
                            selfcontained = T)
    
  }
  
}))


