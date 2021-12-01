library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/STAR/original/coverage/samtools_depth")

# features
features <- import.bed(con = "../L1_6000bp_ordered.bed")

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

plotScatterAndRanges <- function(name_index, x){
  
  bam2binsReduced <- function(bam_subset){
    bam_subset <- bam_subset[which(!width(bam_subset) > 300), ]
    ir <- reduce(ranges(bam_subset))
    bins <- disjointBins(IRanges(start(ir), end(ir) + 1))
    dat <- cbind(as.data.frame(ir), bin = bins)
    if (nrow(dat) == 0){
      dat <- data.frame(start = 0, end = 0, width = 0, bin = 0)
    }
    return(dat)
  }
  
  # subseting coverage for plot
  feature_range <- start(features[name_index, ]) : end(features[name_index, ])
  coverage_plot <- coverage_all_std[coverage_all_std$name == as.character(seqnames(features[name_index, ])) & 
                                    coverage_all_std$pos %in% feature_range, ]
  coverage_plot <- gather(coverage_plot[, 2:5], sample, count, -pos)
  coverage_plot <- coverage_plot[!coverage_plot$count == 0, ]
  
  # subseting .bam for plot
  bam_x4_subset <- granges(subsetByOverlaps(bam_x4, features[name_index, ], minoverlap = 100))
  bam_x5_subset <- granges(subsetByOverlaps(bam_x5, features[name_index, ], minoverlap = 100))
  bam_x6_subset <- granges(subsetByOverlaps(bam_x6, features[name_index, ], minoverlap = 100))

  bam_x4_bins_reduced <- bam2binsReduced(bam_x4_subset)
  bam_x4_bins_reduced$sample <- "11919X4"
  bam_x4_bins_reduced$bin <- -x 
  
  bam_x5_bins_reduced <- bam2binsReduced(bam_x5_subset)
  bam_x5_bins_reduced$sample <- "11919X5"
  bam_x5_bins_reduced$bin <- -x 
  
  bam_x6_bins_reduced <- bam2binsReduced(bam_x6_subset)
  bam_x6_bins_reduced$sample <- "11919X6"
  bam_x6_bins_reduced$bin <- -x 
  bam_bins_reduced <- rbind(bam_x4_bins_reduced, bam_x5_bins_reduced, bam_x6_bins_reduced)

  current_plot <- ggplot() + 
    geom_point(data = coverage_plot, aes(x = pos, y = count, color = sample), size = 0.3) +
    geom_rect(data = bam_bins_reduced, color = "black", aes(xmin = start, xmax = end, ymin = bin, ymax = bin + x*0.6, fill = sample)) +
    scale_x_continuous(limits = c(start(features[name_index, ]) - 100, end(features[name_index, ]) + 100),
                       breaks = seq(from = start(features[name_index, ]), to = end(features[name_index, ]), by = 100), 
                       labels = c(paste0("start: ", start(features[name_index, ])),
                                  rep("", length(seq(from = start(features[name_index, ]), to = end(features[name_index, ]),by = 100)) - 2),
                                  paste0("end: ", end(features[name_index, ])))) +
    scale_y_continuous(limits = c(min(bam_bins_reduced$bin), max(coverage_plot$count))) + 
    facet_grid(sample ~ .) +
    ggtitle(paste(mcols(features[name_index, 1])$name, features[name_index, 1], sep = " "))
  return(current_plot)
}

plot1 <- plotScatterAndRanges(name_index = 1, x = 500)
plot2 <- plotScatterAndRanges(name_index = 2, x = 20)
plot3 <- plotScatterAndRanges(name_index = 3, x = 50)
plot4 <- plotScatterAndRanges(name_index = 4, x = 5)
plot5 <- plotScatterAndRanges(name_index = 5, x = 3)
plot_1_5 <- list(plot1, plot2, plot3, plot4, plot5)
plot_all <- lapply(6:8, FUN = plotScatterAndRanges, x = 3)
plot_all <- c(plot_1_5, plot_all)

pdf("std_coverage_scatter_original_mapping.pdf", width = 19.9, height = 10)
invisible(lapply(plot_all, print))
dev.off()

