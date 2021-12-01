library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/STAR/L1_over_6000bp/mismatch_0.01/coverage/filtered/")

# features
features <- import.bed(con = "L1_6000bp_top10.bed")
# start(features) <- start(features) + 200
# end(features) <- end(features) - 200

# coverage samtoolsdepth
coverage_11919X4_std <- read.delim("11919X4_coverage_std.txt", header = F, col.names = c("name", "pos", "count"))
coverage_11919X5_std <- read.delim("11919X5_coverage_std.txt", header = F, col.names = c("name", "pos", "count"))
coverage_11919X6_std <- read.delim("11919X6_coverage_std.txt", header = F, col.names = c("name", "pos", "count"))
coverage_all_std <- full_join(coverage_11919X4_std, coverage_11919X5_std, by = c("name", "pos"))
coverage_all_std <- full_join(coverage_all_std, coverage_11919X6_std, by = c("name", "pos"))
colnames(coverage_all_std) <- c("name", "pos", paste0("11919X", 4:6))
coverage_all_std[is.na(coverage_all_std)] <- 0

# bam files for bin coverage
filenames <- file.path(paste0("../../11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))
bam_x4 <- readGAlignments(filenames[1])
bam_x5 <- readGAlignments(filenames[2])
bam_x6 <- readGAlignments(filenames[3])

plotScatterAndRanges <- function(name_index){
  
  bam2binsReduced <- function(bam_subset){
    bam_subset <- bam_subset[which(!width(bam_subset) > 300), ]
    ir <- reduce(ranges(bam_subset))
    bins <- disjointBins(IRanges(start(ir), end(ir) + 1))
    dat <- cbind(as.data.frame(ir), bin = bins) 
    return(dat)
  }
  
  # subseting coverage for plot
  coverage_plot <- coverage_all_std[coverage_all_std$name == as.character(seqnames(features[name_index, ])), ]
  coverage_plot <- gather(coverage_plot[, 2:5], sample, count, -pos)
  coverage_plot <- coverage_plot[!coverage_plot$count == 0, ]
  
  # subseting .bam for plot
  bam_x4_subset <- granges(subsetByOverlaps(bam_x4, features[name_index, ]))
  bam_x4_subset <- bam_x4_subset[which(!end(bam_x4_subset) > end(features[name_index, ])), ]
  bam_x5_subset <- granges(subsetByOverlaps(bam_x5, features[name_index, ]))
  bam_x5_subset <- bam_x5_subset[which(!end(bam_x5_subset) > end(features[name_index, ])), ]
  bam_x6_subset <- granges(subsetByOverlaps(bam_x6, features[name_index, ]))
  bam_x6_subset <- bam_x6_subset[which(!end(bam_x6_subset) > end(features[name_index, ])), ]

  bam_x4_bins_reduced <- bam2binsReduced(bam_x4_subset)
  bam_x4_bins_reduced$sample <- "11919X4"
  bam_x4_bins_reduced$bin <- -9
  
  bam_x5_bins_reduced <- bam2binsReduced(bam_x5_subset)
  bam_x5_bins_reduced$sample <- "11919X5"
  bam_x5_bins_reduced$bin <- -6
  
  bam_x6_bins_reduced <- bam2binsReduced(bam_x6_subset)
  bam_x6_bins_reduced$sample <- "11919X6"
  bam_x6_bins_reduced$bin <- -3
  
  bam_bins_reduced <- rbind(bam_x4_bins_reduced, bam_x5_bins_reduced, bam_x6_bins_reduced)

  current_plot <- ggplot() + 
    geom_point(data = coverage_plot, aes(x = pos, y = count, color = sample), size = 0.1) +
    # geom_rect(data = bam_bins_reduced, aes(xmin = start, xmax = end, ymin = bin, ymax = bin + 2.9, fill = sample)) +
    scale_x_continuous(limits = c(0, 8000)) + 
    # scale_y_continuous(limits = c(-500, 8500)) +
    ggtitle(features[name_index, 1]) 
  ggsave(filename = paste0(name_index, "_std_coverage_position_not_scaled.png"), plot = current_plot)
  
  return(name_index)
}
sapply(1:10, plotScatterAndRanges)  

