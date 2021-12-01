### INFO: plots interactive MA plot using plotly and htmlwidgets
### DATE: Tue Mar 06 19:39:21 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
# Sys.setenv("plotly_username" = "fihorvat")
# Sys.setenv("plotly_api_key" = "K7MsyxHYdSADpyy9oeIT")

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/Svoboda/mESC_oocytes_2018/results")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(htmlwidgets)
library(plotly)

######################################################## PATH VARIABLES
outpath <- getwd()

######################################################## SOURCE FILES
# lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
# source(file.path(lib_path, "wideScreen.R"))
# source(file.path(lib_path, "headt.R"))
# source(file.path(lib_path, "asdf.R"))
# wideScreen()

######################################################## FUNCTIONS

######################################################## READ DATA
# get results
results_df <- readr::read_csv(file = "diffExp.MII.SOMvsHET.GRCm38.89.all.csv")

######################################################## MAIN CODE
### MA plot
# data for plot
results_plot <- 
  results_df %>% 
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, 
                mean_FPKM_HET, mean_FPKM_SOM) %>% 
  dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                sign = ifelse(padj < 0.1, "yes", "no"),
                regulation = ifelse(lfc > 0, "up", "down"), 
                regulation = replace(regulation, sign == "no", "not_sign"), 
                regulation = factor(regulation, levels = c("not_sign", "up", "down")), 
                gene_description = str_replace(gene_description, " \\[.*", ""), 
                gene_description = replace(gene_description, is.na(gene_description), "")) %>%  
  # dplyr::filter(sign == "yes") %>%
  dplyr::arrange(regulation)

# make interactive MA plot
p <- 
  plot_ly(data = results_plot, 
          x = ~mean, 
          y = ~lfc, 
          text = ~paste("</br> log2FC: ", round(lfc, 3), 
                        "</br> mean exp.: ", round(mean, 3), 
                        "</br> HET avg. FPKM:", mean_FPKM_HET,
                        "</br> SOM avg. FPKM:", mean_FPKM_SOM,
                        "</br>", gene_id,
                        "</br>", gene_name,
                        "</br>", gene_description),
          color = ~regulation, 
          colors = c("gray32", "blue3", "red3"), 
          alpha = 0.75, 
          hoverinfo = "text") %>%
  add_markers() %>%
  layout(xaxis = list(title = "mean expression", type = "log"),
         yaxis = list(title = "log2FoldChange SOM/HET"))

# save to .html
htmlwidgets::saveWidget(as_widget(p), file = file.path(outpath, "MAplot.MII.SOMvsHET.DESeq2.html"), selfcontained = T)


# ### FPKM MA-plot
# fpkm_MA_plot_df <-
#   results_plot %>%
#   dplyr::mutate(log2FC = log2(mean_FPKM_SOM) - log2(mean_FPKM_HET),
#                 mean_exp = (mean_FPKM_SOM + mean_FPKM_HET) / 2)
# 
# # make interactive MA plot
# p2 <-
#   plot_ly(data = fpkm_MA_plot_df,
#           x = ~mean_exp,
#           y = ~log2FC,
#           text = ~paste("</br> log2FC: ", round(log2FC, 3),
#                         "</br> mean exp.: ", round(mean_exp, 3),
#                         "</br>", gene_id,
#                         "</br>", gene_name,
#                         "</br>", gene_description),
#           color = ~regulation,
#           colors = c("gray32", "blue3", "red3"),
#           alpha = 0.75,
#           hoverinfo = "text") %>%
#   add_markers() %>%
#   layout(xaxis = list(title = "mean expression", type = "log"),
#          yaxis = list(title = "log2FoldChange SOM/HET"))
# 
# # save to .html
# htmlwidgets::saveWidget(as_widget(p2), file = file.path(outpath, "MAplot.MII.SOMvsHET.FPKM.html"), selfcontained = T)

