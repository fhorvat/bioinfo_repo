library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/STAR/L1_over_6000bp/mismatch_0.01/coverage/filtered/samtools_depth_raw")
bed_top10 <- read.delim("../L1_6000bp_top10.bed", header = F, stringsAsFactors = F)

# coverage samtoolsdepth
coverage_11919X4_std <- read.delim("11919X4_coverage_std.txt", header = F, col.names = c("name", "pos", "count"))
coverage_11919X5_std <- read.delim("11919X5_coverage_std.txt", header = F, col.names = c("name", "pos", "count"))
coverage_11919X6_std <- read.delim("11919X6_coverage_std.txt", header = F, col.names = c("name", "pos", "count"))
coverage_all_std <- full_join(coverage_11919X4_std, coverage_11919X5_std, by = c("name", "pos"))
coverage_all_std <- full_join(coverage_all_std, coverage_11919X6_std, by = c("name", "pos"))
colnames(coverage_all_std) <- c("name", "pos", paste0("11919X", 4:6))
coverage_all_std[is.na(coverage_all_std)] <- 0

# plot
coverage_plot <- function(name_index){
  coverage_plot <- coverage_all_std[coverage_all_std$name == bed_top10[name_index, 1], ]
  coverage_plot <- gather(coverage_plot[, 2:5], sample, count, -pos)
  current_plot <- ggplot(coverage_plot, aes(x = pos, y = count, color = sample)) + 
    geom_point(size = 0.1) +
    scale_x_continuous(limit = c(175, bed_top10[name_index, 3] - 175)) +
    ggtitle(paste0(bed_top10[name_index, 1], "_std"))
  ggsave(filename = paste0(name_index, "_std.png"), plot = current_plot)
  return(name_index)
}

sapply(1:10, coverage_plot)
