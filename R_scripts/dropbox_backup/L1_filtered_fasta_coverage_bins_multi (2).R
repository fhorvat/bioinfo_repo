library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(DESeq2)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/")
features <- import.bed(con = "STAR/filtered_fasta/Pepa_filter/CNOT6_WT_GV/coverage/full_length_L1_edit_ordered.bed")
start(features) <- 1

# bam_files
filenames <- file.path(paste0("STAR/filtered_fasta/Pepa_filter/CNOT6_WT_GV/11919X", 4:6, "_Aligned.sortedByCoord.out.bam"))
bam_x4 <- readGAlignments(filenames[1])
bam_x5 <- readGAlignments(filenames[2])
bam_x6 <- readGAlignments(filenames[3])

plotRangesMulti <- function(name_index){

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
  
  bam_x4_subset <- granges(subsetByOverlaps(bam_x4, features[name_index, ]))
  bam_x5_subset <- granges(subsetByOverlaps(bam_x5, features[name_index, ]))
  bam_x6_subset <- granges(subsetByOverlaps(bam_x6, features[name_index, ]))
  
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

  x <- max(bam_bins$bin) / 20
  if (x < 2){
    x <- 2
  }
  
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
    scale_x_continuous(limits = c(-101, end(features[name_index, ]) + 101), 
                       breaks = seq(0, to = end(features[name_index, ]), by = 100)) +
    scale_y_continuous(limits = c(min(bam_bins_reduced$bin), max(bam_bins$bin) + 3)) + 
    facet_grid(sample ~ .) +
    ggtitle(seqnames(features[name_index, ])) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(current_plot)
}

# plot1 <- plotRangesMulti(name_index = 1, x = 5000)
# plot2 <- plotRangesMulti(name_index = 2, x = 2000)
# plot3 <- plotRangesMulti(name_index = 3, x = 1000)
# plot4 <- plotRangesMulti(name_index = 4, x = 500)
# plot5 <- plotRangesMulti(name_index = 5, x = 200)
# plot6 <- plotRangesMulti(name_index = 6, x = 100)

plot_all <- lapply(1:length(features), FUN = plotRangesMulti)

pdf("STAR/filtered_fasta/Pepa_filter/full_length_L1_coverage_bins_CNOT6.pdf", width = 19.9, height = 10)
invisible(lapply(plot_all, print))
dev.off()


