library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/STAR/filtered_fasta/mismatch_0.01/coverage")
features <- import.bed(con = "L1_6000bp_filtered_fasta_ordered.bed")
  
# removing 200bp flanking regions from bed features
start(features) <- start(features) + 199
end(features) <- end(features) - 199

# bam files for bin coverage
filenames <- file.path(paste0("../11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))
bam_x4 <- readGAlignments(filenames[1])
bam_x5 <- readGAlignments(filenames[2])
bam_x6 <- readGAlignments(filenames[3])

plotRangesMulti <- function(name_index, x){

  bam2bins <- function(bam_subset){
    bam_subset <- bam_subset[which(!width(bam_subset) > 300), ]
    ir <- ranges(bam_subset)
    bins <- disjointBins(IRanges(start(ir), end(ir) + 1)) 
    dat <- cbind(as.data.frame(ir), bin = bins)
    if (nrow(dat) == 0){
      dat <- data.frame(start = 0, end = 0, width = 0, bin = 0)
    }
    return(dat)
  }
  
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
  
  bam_x4_subset <- granges(subsetByOverlaps(bam_x4, features[name_index, ], minoverlap = 100))
  bam_x5_subset <- granges(subsetByOverlaps(bam_x5, features[name_index, ], minoverlap = 100))
  bam_x6_subset <- granges(subsetByOverlaps(bam_x6, features[name_index, ], minoverlap = 100))
  bam_x4_bins <- bam2bins(bam_x4_subset)
  bam_x4_bins$sample <- "11919X4"
  bam_x4_bins$bin <- seq(0, nrow(bam_x4_bins), 1)[1:nrow(bam_x4_bins)] 
  bam_x5_bins <- bam2bins(bam_x5_subset)
  bam_x5_bins$sample <- "11919X5"
  bam_x5_bins$bin <- seq(0, nrow(bam_x5_bins), 1)[1:nrow(bam_x5_bins)] 
  bam_x6_bins <- bam2bins(bam_x6_subset)
  bam_x6_bins$sample <- "11919X6"
  bam_x6_bins$bin <- seq(0, nrow(bam_x6_bins), 1)[1:nrow(bam_x6_bins)] 
  
  bam_bins <- rbind(bam_x4_bins, bam_x5_bins, bam_x6_bins)

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
    geom_rect(data = bam_bins, aes(xmin = start, xmax = end, ymin = bin, ymax = bin + 0.8, fill = sample)) +
    geom_rect(data = bam_bins_reduced, color = "black", aes(xmin = start, xmax = end, ymin = bin, ymax = bin + x*0.8, fill = sample)) +
    scale_x_continuous(limits = c(-100, end(features[name_index, ])), 
                       breaks = seq(from = start(features[name_index, ]), to = end(features[name_index, ]) + 200, by = 100),
                       labels = seq(0, end(features[name_index, ]), by = 100)) +
    scale_y_continuous(limits = c(min(bam_bins_reduced$bin), max(bam_bins$bin) + 3)) + 
    facet_grid(sample ~ .) +
    ggtitle(paste(mcols(features[name_index, ])$name, seqnames(features[name_index, ]), sep = " ")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(current_plot)
}

plot1 <- plotRangesMulti(name_index = 1, x = 2000)
plot2 <- plotRangesMulti(name_index = 2, x = 500)
plot3 <- plotRangesMulti(name_index = 3, x = 200)
plot4 <- plotRangesMulti(name_index = 4, x = 50)
plot5 <- plotRangesMulti(name_index = 5, x = 10)
plot6 <- plotRangesMulti(name_index = 6, x = 2)
plot_1_6 <- list(plot1, plot2, plot3, plot4, plot5, plot6)

plot_all <- lapply(7:1000, FUN = plotRangesMulti, x = 2)
plot_all <- c(plot_1_6, plot_all)

pdf("L1_6000bp_filtered_fasta_mapping_coverage_bins_1000.pdf", width = 19.9, height = 10)
invisible(lapply(plot_all, print))
dev.off()


