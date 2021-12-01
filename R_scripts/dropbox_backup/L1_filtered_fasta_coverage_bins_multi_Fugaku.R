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
features <- import.bed(con = "STAR/filtered_fasta/Pepa_filter/full_length_L1_Fugaku_testis_edit_ordered.bed")
start(features) <- 1

# bam_files
filenames <- file.path(c("STAR/filtered_fasta/Pepa_filter/Fugaku_blastocyst/s_Blast.WE_Aligned.sortedByCoord.out.bam", 
                         "STAR/filtered_fasta/Pepa_filter/Fugaku_GV/s_GV.WE_Aligned.sortedByCoord.out.bam", 
                         "STAR/filtered_fasta/Pepa_filter/Encode_testis/testis_Aligned.sortedByCoord.out.bam"))

bam_Blast <- readGAlignments(filenames[1])
bam_GV <- readGAlignments(filenames[2])
bam_testis <- readGAlignments(filenames[3])

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
  
  bam_Blast_subset <- granges(subsetByOverlaps(bam_Blast, features[name_index, ]))
  bam_GV_subset <- granges(subsetByOverlaps(bam_GV, features[name_index, ]))
  bam_testis_subset <- granges(subsetByOverlaps(bam_testis, features[name_index, ]))
  
  bam_Blast_bins <- bam2bins(bam_Blast_subset)
  bam_Blast_bins$sample <- "Fugaku_blastocyst"
  bam_Blast_bins$bin <- seq(0, nrow(bam_Blast_bins), 1)[1:nrow(bam_Blast_bins)] 
  
  bam_GV_bins <- bam2bins(bam_GV_subset)
  bam_GV_bins$sample <- "Fugaku_GV"
  bam_GV_bins$bin <- seq(0, nrow(bam_GV_bins), 1)[1:nrow(bam_GV_bins)] 
  
  bam_testis_bins <- bam2bins(bam_testis_subset)
  bam_testis_bins$sample <- "Encode_testis"
  bam_testis_bins$bin <- seq(0, nrow(bam_testis_bins), 1)[1:nrow(bam_testis_bins)] 
  
  bam_bins <- rbind(bam_Blast_bins, bam_GV_bins, bam_testis_bins)

  x <- max(bam_bins$bin) / 20
  if (x < 2){
    x <- 2
  }
  
  bam_Blast_bins_reduced <- bam2binsReduced(bam_Blast_subset)
  bam_Blast_bins_reduced$sample <- "Fugaku_blastocyst"
  bam_Blast_bins_reduced$bin <- -x 
  bam_GV_bins_reduced <- bam2binsReduced(bam_GV_subset)
  bam_GV_bins_reduced$sample <- "Fugaku_GV"
  bam_GV_bins_reduced$bin <- -x 
  bam_testis_bins_reduced <- bam2binsReduced(bam_testis_subset)
  bam_testis_bins_reduced$sample <- "Encode_testis"
  bam_testis_bins_reduced$bin <- -x 
  
  bam_bins_reduced <- rbind(bam_Blast_bins_reduced, bam_GV_bins_reduced, bam_testis_bins_reduced)
  
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

plot1 <- plotRangesMulti(name_index = 1)
plot2 <- plotRangesMulti(name_index = 2)
plot3 <- plotRangesMulti(name_index = 3)
plot4 <- plotRangesMulti(name_index = 4)
plot5 <- plotRangesMulti(name_index = 5)
plot6 <- plotRangesMulti(name_index = 6)

plot_all <- lapply(1:length(features), FUN = plotRangesMulti)

pdf("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/STAR/filtered_fasta/Pepa_filter/full_length_L1_coverage_bins_Fugaku_testis.pdf", width = 19.9, height = 10)
invisible(lapply(plot_all, print))
dev.off()
