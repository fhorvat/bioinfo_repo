### INFO: 
### DATE: Tue Feb 04 18:21:56 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/Analysis/DicerX_Tarbp_comparison")

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

library(plotly)
library(htmlwidgets)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# DicerX FPM path
dicerx_fpm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Analysis/expression/miRBase/miRBase.22.mm10.20181605.FPM_mean.csv"

# Tarbp2 FPM path
tarbp2_fpm_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Pullagura_2018_Genetics_PRJNA423238/Analysis/expression/miRBase/miRBase.22.mm10.20181605.FPM_mean.csv"

######################################################## READ DATA
# read DicerX FPM
dicerx_fpm <- readr::read_csv(dicerx_fpm_path) 

# read Tarbp2 FPM
tarbp2_fpm <- readr::read_csv(tarbp2_fpm_path)

######################################################## MAIN CODE
# clean files, join
fpm_tb <- 
  dicerx_fpm %>%
  dplyr::select(gene_id, DicerX_KO = KO, DicerX_WT = WT, coordinates) %>%
  dplyr::left_join(., tarbp2_fpm %>% dplyr::select(gene_id, Tarbp2_Mut, Tarbp2_WT), by = "gene_id") %>% 
  dplyr::mutate(DicerX_log2_ratio = log2((DicerX_KO + 0.1) / (DicerX_WT + 0.1)), 
                Tarbp2_log2_ratio = log2((Tarbp2_Mut + 0.1) / (Tarbp2_WT + 0.1))) %>%   
  dplyr::select(gene_id, DicerX_log2_ratio, Tarbp2_log2_ratio, DicerX_KO, DicerX_WT, Tarbp2_Mut, Tarbp2_WT, coordinates)

# save
fpm_tb %>% 
  dplyr::select(gene_id, coordinates, 
                DicerX_KO_FPM = DicerX_KO, DicerX_WT_FPM = DicerX_WT, Tarbp2_Mut_FPM = Tarbp2_Mut, Tarbp2_WT_FPM = Tarbp2_WT, 
                DicerX_log2_ratio, Tarbp2_log2_ratio) %T>% 
  readr::write_csv(., file.path(outpath, str_c("miRBase.22.mm10.20181605.FPM_mean", "DicerX_vs_Tarbp2.log2_ratio.FPM", "csv", sep = ".")))

# ### compare two wild-types
# # build a linear model
# linearMod <- lm(DicerX_WT ~ Tarbp2_WT, data = fpm_tb)  # build linear regression model on full data
# summary(linearMod)
# 
# # plot 
# ln_plot <- 
#   ggplot() +
#   geom_point(data = fpm_tb, 
#              mapping = aes(x = DicerX_WT, y = Tarbp2_WT), alpha = 1, size = 3) +
#   geom_smooth(data = fpm_tb,
#               mapping = aes(x = DicerX_WT, y = Tarbp2_WT),
#               method = lm) +
#   scale_x_continuous(limits = c(0, 500)) +
#   scale_y_continuous(limits = c(0, 500)) +
#   xlab("DicerX WT FPM") +
#   ylab("Tarbp2 WT FPM") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         legend.position = "none")
# 
# # save
# ggsave(filename = file.path(outpath, str_c("miRBase.22.mm10.20181605.FPM_mean", "DicerX_vs_Tarbp2_WT.linear_regression", "png", sep = ".")), 
#        plot = ln_plot, 
#        height = 10, width = 10)

# get limits
axis_limits <- 
  c(fpm_tb$DicerX_log2_ratio, fpm_tb$Tarbp2_log2_ratio) %>% 
  abs(.) %>% 
  max(.) %>% 
  ceiling(.)

# crosshair plot
cross_plot <-
  ggplot(fpm_tb, aes(x = DicerX_log2_ratio, y = Tarbp2_log2_ratio)) +
  geom_point(shape = 16, size = 3, color = "gray30", alpha = 0.75) +
  scale_x_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  scale_y_continuous(limits = c(-axis_limits, axis_limits), breaks = seq(-axis_limits, axis_limits, 2)) +
  # scale_color_manual(values = c(lnc1 = "red3", dormant = "#1a75ff", histone = "black", protein_coding = "black")) +
  # scale_size_manual(values = c(lnc1 = 4, dormant = 2.5, histone = 2.5, protein_coding = 2.5)) +
  # scale_alpha_manual(values = c(lnc1 = 1, dormant = 1, histone = 0.2, protein_coding = 0.2)) + 
  # guides(color = FALSE, size = FALSE, alpha = FALSE) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab(str_c("log2 (DicerX KO / WT) FPM")) +
  ylab(str_c("log2 (Tarbp2 KO / WT) FPM")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = 0.3),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("miRBase.22.mm10.20181605.FPM_mean", "DicerX_vs_Tarbp2.log2_ratio.crosshair", "png", sep = ".")),
       plot = cross_plot,
       width = 10, height = 10)


### interactive
# interactive plot
interactive_plot <-
  plotly::plot_ly(data = fpm_tb,
                  x = ~DicerX_log2_ratio,
                  y = ~Tarbp2_log2_ratio,
                  text = ~paste("</br>", gene_id,
                                "</br> DicerX KO FPM: ", round(DicerX_KO, 2),
                                "</br> DicerX WT FPM: ", round(DicerX_WT, 2),
                                "</br> Tarbp2 Mut FPM: ", round(Tarbp2_Mut, 2),
                                "</br> Tarbp2 WT FPM: ", round(Tarbp2_WT, 2),
                                "</br>", coordinates),
                  alpha = 0.75,
                  hoverinfo = "text") %>%
  add_markers() %>%
  layout(xaxis = list(title = "log2 (DicerX KO / WT) FPM"),
         yaxis = list(title = "log2 (Tarbp2 KO / WT) FPM"))

# save as html widget
htmlwidgets::saveWidget(plotly::as_widget(interactive_plot),
                        file = file.path(outpath, str_c("miRBase.22.mm10.20181605.FPM_mean", "DicerX_vs_Tarbp2.log2_ratio.crosshair", "html", sep = ".")), 
                        selfcontained = T)

