### INFO: summarizes summarizedExperiments, calculates RPM
### DATE: Mon Feb 25 14:20:36 2019
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
library(openxlsx)
library(purrr)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)

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

# sample table path
sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/developmental_stages.RNAseq.sample_table.20190429.csv"

# datasets path
datasets_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/developmental_stages/datasets"

# experiment list
experiment_list <- c("Deng_2014_Science_GSE45719", 
                     # "Smallwood_2011_NatGenet_PRJEB2547",
                     # "Hamazaki_2015_Development_PRJDB2994", 
                     "Veselovska_2015_GenomeBiol_GSE70116",
                     "Yamaguchi_2013_CellRes_GSE41908", 
                     "Gan_2013_NatCommun_GSE35005", 
                     "ENCODE_2014_Nature_GSE49417")

# perfect read stats path
read_stats_paths <- file.path(datasets_path, experiment_list, "count_perfect_reads")

# read stats list 
read_stats_list <- list.files(read_stats_paths, pattern = ".*\\.read_stats.txt", full.names = T)

# counts path 
read_counts_paths <- file.path(datasets_path, experiment_list, "read_counts")

# read counts list 
read_counts_list <- list.files(read_counts_paths, pattern = "rmsk_counts\\..*\\.RDS", full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- 
  read_csv(sample_table_path) %>% 
  dplyr::filter(experiment %in% experiment_list) %>% 
  dplyr::select("sample_id", "stage") %>% 
  mutate(sample_id = str_remove(sample_id, "\\.total\\.bam")) %>% 
  as.data.table(.)

# read and sum read stats
read_stats <- 
  map(read_stats_list, function(path){
    
    # read and clean
    read_delim(path, delim = "\t") %>% 
      dplyr::mutate(genome.mapped_minus_rDNA = ifelse(str_detect(sample_id, "\\.PE"),
                                                      round(genome.mapped_minus_rDNA / 2),
                                                      genome.mapped_minus_rDNA)) %>%
      dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>%
      mutate(sample_id = str_remove(sample_id, "\\.total\\.bam")) %>% 
      mutate(library_size = round(library_size / 1E6, 6))
    
  }) %>% 
  bind_rows(.) %>% 
  as.data.table(.)

# read counts
read_counts <- 
  map(read_counts_list, readRDS) %>% 
  bind_rows(.) %>% 
  as.data.table(.)

######################################################## MAIN CODE
### clean stage names and repeats
# clean stage names
stage_clean <- 
  tibble(stage = c("PGC_9.5", "PGC_11.5", "PGC_13.5_m", "PGC_13.5_f", 
                   "primitive_SG_A", "SG_B", "leptotene_SC", "pachytene_SC", "round_ST", 
                   "nonGrowing_oocytes", "growing_oocytes_d8_14", "GV_oocytes",
                   "zygote", "2C", "4C", "8C", "16C", "mid_blast",
                   "placenta_adult8wks"), 
         stage_clean = c("PGC 9.5", "PGC 11.5", "male PGC 13.5", "female PGC 13.5", 
                         "primitive SG A", "SG B", "leptotene SC", "pachytene SC", "round ST",
                         "non-growing oocyte", "growing oocyte days 8-14", "GV ",
                         "zygote", "2-cell", "4-cell", "8-cell", "16-cell", "blastocyst 94h", 
                         "placenta")) %>%
  mutate(stage_clean = factor(stage_clean, levels = stage_clean))

# clean repeatMasker classes
rmsk_clean <- 
  tibble(rmsk_id = c("LINE1_all", "L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F", 
                     "IAP_all", "IAPEy_int_IAPLTR1_Mm", "MERVL_int_MT2_Mm", "MTA_int_LTR", "ORR1A_int_LTR"), 
         rmsk_clean = c("LINE 1 all", "L1Md_T", "L1Md_A", "L1Md_Gf", "L1Md_F2", "L1_Mus2", "L1_Mus1", "L1Md_F", 
                        "IAP all", "IAPEy_int + IAPLTR1_Mm", "MERVL_int + MT2_Mm", "MTA_int + LTRs", "ORR1A_int + LTRs")) %>% 
  mutate(rmsk_clean = factor(rmsk_clean, levels = rmsk_clean))


### calculate RPM, plot as heatmap
# get average RPM
avg_rpm_tb <- 
  copy(read_counts) %>% 
  .[read_stats, on = "sample_id", `:=`(rpm = round((count / i.library_size), 3))] %>% 
  .[sample_table, on = "sample_id", `:=`(stage = i.stage)] %>% 
  .[, list(avg_rpm = mean(rpm)), by = .(rmsk_id, stage)] %>% 
  .[] %>% 
  as.tibble(.)

# get table for heatmap
rpm_heatmap_tb <- 
  avg_rpm_tb %>% 
  mutate(log10_rpm = log10(avg_rpm)) %>%
  right_join(., stage_clean, by = "stage") %>% 
  right_join(., rmsk_clean, by = "rmsk_id") %>% 
  dplyr::select(stage_clean, rmsk_clean, log10_rpm) %>%
  tidyr::spread(key = stage_clean, value = log10_rpm) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "rmsk_clean") %>% 
  as.matrix(.)

# get columns and rows to put gaps
which.gap.col <- which(colnames(rpm_heatmap_tb) %in% c("female PGC 13.5", "round ST", "fully-grown oocyte", "GV ", "blastocyst 94h", "placenta"))
which.gap.row <- which(rownames(rpm_heatmap_tb) %in% c("L1Md_F", "IAPEy_int + IAPLTR1_Mm", "MERVL_int + MT2_Mm", "MTA_int + LTRs", "ORR1A_int + LTRs"))

# set colors and breaks
my.cols <- c("white", "white", "white", "white", str_c("grey", 99:1), "black", "black")
my.breaks <- c(-1, 0, seq(1, 5.5, length = 100), 6)

# my.cols <- c("white", str_c("grey", 99:1), "black")
# my.breaks <- seq(0, 5, length = 101)

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
                   filename = file.path(outpath, str_c("heatmap.rmsk.L1_and_LTRs.20190429.RPM.png")),
                   height = 20,
                   width = 20)

