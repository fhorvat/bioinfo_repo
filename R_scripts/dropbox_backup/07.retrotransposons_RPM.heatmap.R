### INFO: make heatmap from retrotransposon expression
### DATE: Tue Feb 26 23:34:26 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages")

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

library(pheatmap)
library(RColorBrewer)
library(ggthemes)
library(grid)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# changes draw_legend function in pheatmap
my.draw.legend <-  function (color, breaks, legend, ...){
  
  library(grid)
  
  height = min(unit(1, "npc"), unit(150, "bigpts"))
  legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos = height * legend_pos + (unit(1, "npc") - height)
  breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
  breaks = height * breaks + (unit(1, "npc") - height)
  h = breaks[-1] - breaks[-length(breaks)]
  h2 = breaks[length(breaks)] - breaks[1]
  rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = NA))
  rect1 = rectGrob(x = 0, y = breaks[1], width = unit(10, "bigpts"), height = h2, hjust = 0, vjust = 0, gp = gpar(fill = NA, col = "grey60"))
  text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
  res = grobTree(rect, text, rect1)
  
  return(res)
  
}

# change pheatmap:::draw_legendwith my.draw.legend
assignInNamespace("draw_legend", my.draw.legend, ns = "pheatmap")

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# summarizedExperiments path
se_path <- file.path(inpath, "summarizedExperiments.all_LINEs_LTRs")

######################################################## READ DATA
# experiment list
experiment_list <- c(
  "Hamazaki_2015_Development_PRJDB2994", 
  "Smallwood_2011_NatGenet_PRJEB2547", 
  "Yamaguchi_2013_CellRes_GSE41908",
  "Gan_2013_NatCommun_GSE35005", 
  "Deng_2014_Science_GSE45719",
  "ENCODE_2014_Nature_GSE49417"
)

# list all RPM retrotransposon table
retro_rpm_list <- list.files(se_path, pattern = "rmsk\\.L1_and_LTRs\\..*\\.avgRPM\\.RDS", full.names = T)

# read RPM 
retro_rpm <- 
  map(retro_rpm_list, readRDS) %>% 
  set_names(., str_extract(retro_rpm_list, str_c(experiment_list, collapse = "|")))

######################################################## MAIN CODE
# transpose RPMs list, join to one data.table per filter
retro_rpm <- 
  lapply(c("whole", "removed_exons", "removed_genes"), function(filter_name){
    
    # get experiment table
    lapply(retro_rpm, "[[", filter_name) %>% 
      do.call(rbind, .)
    
  }) %>% 
  set_names(., c("whole", "removed_exons", "removed_genes"))

# clean stage names
stage_clean <- 
  tibble(stage = c("PGC_9.5", "PGC_11.5", "PGC_13.5_m",
                   "primitive_SG_A", "SG_A", "SG_B", "leptotene_SC", "pachytene_SC", "round_ST", "elongative_ST", 
                   "PGC_13.5_f", "d10_oocyte", "FG_GV_oocyte", "MII", "zygote", "2C", "4C", "8C", "16C", 
                   "early_blast", "mid_blast", "late_blast", "placenta_adult8wks"), 
         stage_clean = c("PGC 9.5", "PGC 11.5", "male PGC 13.5",
                         "primitive SG A", "SG A", "SG B", "leptotene SC", "pachytene SC", "round ST", "elongative ST", 
                         "female PGC 13.5", "small oocyte", "fully-grown oocyte", "unfertilized egg", 
                         "zygote", "2-cell", "4-cell", "8-cell", "16-cell", 
                         "blastocyst 88h", "blastocyst 94h", "blastocyst 102h", "placenta")) %>%
  mutate(stage_clean = factor(stage_clean, levels = stage_clean))

### clean classes
# LINE1 row like in GenRes
# Line1 subfamilies = L1Md_T L1Md_A L1Md_Gf L1Md_F2 L1_Mus2 L1_Mus1 L1Md_F
# IAP all
# IAPEz-Int with LTRs IAPLTR1_Mm
# MuERV-int with MT2 LTRs = MERVL + MT2_Mm
# MTA all LTRs and ints
# ORR1A all ints and LTR subgroups
retro_clean <- 
  tibble(class = c("LINE1_all", "L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F", 
                   "IAP_all", "IAPEy_int_IAPLTR1_Mm", "MERVL_int_MT2_Mm", "MTA_int_LTR", "ORR1A_int_LTR"), 
         retro_clean = c("LINE 1 all", "L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F", 
                         "IAP all", "IAPEy_int + IAPLTR1_Mm", "MERVL_int + MT2_Mm", "MTA_int + LTRs", "ORR1A_int + LTRs")) %>% 
  mutate(retro_clean = factor(retro_clean, levels = retro_clean))

### plot heatmap for 3 different repeatMasker filters 
# loop through filters
map(names(retro_rpm), function(filter_name){
  
  # read and bind all retrotransposon RPMs
  rpm_heatmap_tb <- 
    retro_rpm[[filter_name]] %>% 
    as.tibble(.) %>% 
    filter(!((stage == "2C") & (experiment == "Hamazaki")),
           stage != "Sertoli") %>% 
    # mutate(log10_rpm = log10(sum_avg_rpm + 1)) %>%
    mutate(log10_rpm = sum_avg_rpm) %>%
    left_join(., stage_clean, by = "stage") %>% 
    left_join(., retro_clean, by = "class") %>% 
    dplyr::select(stage_clean, retro_clean, log10_rpm) %>%
    tidyr::spread(key = stage_clean, value = log10_rpm) %>% 
    as.data.frame(.) %>% 
    tibble::column_to_rownames(., var = "retro_clean") %>% 
    as.matrix(.)
  
  # get columns and rows to put gaps
  which.gap.col <- which(colnames(rpm_heatmap_tb) %in% c("PGC 11.5", "elongative ST", "unfertilized egg", "blastocyst 102h", "placenta"))
  which.gap.row <- which(rownames(rpm_heatmap_tb) %in% c("L1Md_F", "IAPEy_int + IAPLTR1_Mm", "MERVL_int + MT2_Mm", "MTA_int + LTRs", "ORR1A_int + LTRs"))
    
  # set colors and breaks
  my.cols <- c("white", "white", "white", str_c("grey", 99:1), "black", "black")
  my.breaks <- c(0, seq(1.5, 3.5, length = 100), 4)
  
  # plot heatmap with annotation
  pheatmap::pheatmap(rpm_heatmap_tb,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     fontsize_row = 20, 
                     fontsize_col = 20,
                     # col = colorRampPalette(brewer.pal(9, "Greys"))(20),
                     color = my.cols,
                     breaks = my.breaks,
                     cellwidth = 50, 
                     cellheight = 50, 
                     gaps_col = which.gap.col, 
                     gaps_row = which.gap.row,
                     filename = file.path(outpath, str_c("heatmap.rmsk.L1_and_LTRs.20190306", filter_name, "RPM.breaks.png", sep = ".")),
                     height = 20,
                     width = 20)
  
})

