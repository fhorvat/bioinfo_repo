coverageSinglePlot <- function(x, coverage_stage, bam_stage_name){
  MT2SampleGR <- MT2_filtered_ranges_fpkm_filtered_grList_all[[bam_stage_name]][[x]]
  MT2SampleRL <- as(MT2SampleGR, "RangesList")
  MT2SampleCoverage <- Views(coverage_stage, MT2SampleRL)
  MT2SampleCoverage <- MT2SampleCoverage[[as.character(seqnames(MT2SampleGR[1]))]]
  MT2SampleCoverageCounts <- as.vector(unlist(viewApply(MT2SampleCoverage, as.vector)))
  MT2SampleCoveragePos <- unlist(apply(data.frame(x = start(MT2SampleCoverage), y = end(MT2SampleCoverage)), 1, function(x) x[1]:x[2]))
  
  # original sample data
  MT2SampleOriginal <- MT2Top100[MT2Top100$repNames == names(MT2_filtered_ranges_fpkm_filtered_grList_all[[bam_stage_name]][[x]])[1], ]
  MT2SampleExpanded <- MT2Expanded[MT2Expanded$repNames == names(MT2_filtered_ranges_fpkm_filtered_grList_all[[bam_stage_name]][[x]])[1], ]
  
  # counts with names = positions
  names(MT2SampleCoverageCounts) <- MT2SampleCoveragePos
  
  # data.frame for plot
  countsZero <- rep(0, 4e5 + 1)
  names(countsZero) <- MT2SampleExpanded$start : MT2SampleExpanded$end
  MT2SampleCoverageCountsFull <- replace(countsZero, names(countsZero) %in% names(MT2SampleCoverageCounts), MT2SampleCoverageCounts)
  
  MT2SampleCoverageDFPlot <- data.frame(seqnames = seqnames(MT2SampleGR[1]), 
                                        pos = as.numeric(names(MT2SampleCoverageCountsFull)),
                                        counts = MT2SampleCoverageCountsFull,
                                        strand = MT2SampleOriginal$strand,
                                        repName = names(MT2ExpandedFilteredGRList[x]))
  MT2SampleCoverageDFPlot$element <- ifelse((MT2SampleCoverageDFPlot$pos >= MT2SampleOriginal$start & MT2SampleCoverageDFPlot$pos <= MT2SampleOriginal$end), "MT2", "out")
  MT2SampleCoverageDFPlot$pos <- MT2SampleCoverageDFPlot$pos - MT2SampleExpanded$start - 2e5
  if (MT2SampleCoverageDFPlot[1, ]$strand == "-"){
    MT2SampleCoverageDFPlot$pos <- MT2SampleCoverageDFPlot$pos + (MT2SampleOriginal$end - MT2SampleOriginal$start)
    MT2SampleCoverageDFPlot$pos <- rev(MT2SampleCoverageDFPlot$pos)
  }
  MT2SampleCoverageDFPlot$count[MT2SampleCoverageDFPlot$count == 0] <- NA
  MT2SampleCoverageDFPlot$pos[is.na(MT2SampleCoverageDFPlot$count)] <- NA
  
  # plot
  MT2SampleCoveragePlot <- ggplot(MT2SampleCoverageDFPlot, aes(x = pos, y = counts)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = counts, fill = element)) +
    scale_x_continuous(limits = c(-2e5, 2e5)) + 
#     scale_y_continuous(limits = c(0, 3500)) + 
    scale_fill_manual(values = c("red", "black")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_blank()) +
    coord_cartesian(ylim = c(0, 200)) +
    ggtitle(MT2SampleCoverageDFPlot$repName[1]) 
  return(MT2SampleCoveragePlot)
}

allMT2CoverageCounts_GV_plot <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_GV"), 
                                  function(X) coverageSinglePlot(X, coverage_GV, "bam_GV")) 
allMT2CoverageCounts_1C_plot <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_1C"), 
                                  function(X) coverageSinglePlot(X, coverage_1C, "bam_1C")) 
allMT2CoverageCounts_2C_plot <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_2C"), 
                                  function(X) coverageSinglePlot(X, coverage_2C, "bam_2C")) 
allMT2CoverageCounts_2C_aphi_plot <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_2C_aphi"), 
                                       function(X) coverageSinglePlot(X, coverage_2C_aphi, "bam_2C_aphi")) 
allMT2CoverageCounts_4C_plot <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_4C"), 
                                  function(X) coverageSinglePlot(X, coverage_4C, "bam_4C")) 

# plot_list <- invisible(lapply(X = 1:length(MT2ExpandedFilteredGRList), function(X) coverageSinglePlot(X))) 
pdf("MT2_94_GV_fpkm_250_filtered.pdf", width = 10, height = 5)
invisible(lapply(allMT2CoverageCounts_GV_plot, print))
dev.off()

pdf("MT2_94_1C_fpkm_250_filtered.pdf", width = 10, height = 5)
invisible(lapply(allMT2CoverageCounts_1C_plot, print))
dev.off()

pdf("MT2_94_2C_fpkm_250_filtered.pdf", width = 10, height = 5)
invisible(lapply(allMT2CoverageCounts_2C_plot, print))
dev.off()

pdf("MT2_94_2C_aphi_fpkm_250_filtered.pdf", width = 10, height = 5)
invisible(lapply(allMT2CoverageCounts_2C_aphi_plot, print))
dev.off()

pdf("MT2_94_4C_fpkm_250_filtered.pdf", width = 10, height = 5)
invisible(lapply(allMT2CoverageCounts_4C_plot, print))
dev.off()



