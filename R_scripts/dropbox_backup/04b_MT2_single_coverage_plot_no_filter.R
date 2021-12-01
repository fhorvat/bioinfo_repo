MT2_top <- read.csv("MT2_ORR_single_elements_range.csv", stringsAsFactors = F, header = T)
MT2OriginalGR <- makeGRangesFromDataFrame(MT2_top, keep.extra.columns = T)

# expand ranges 200kb up and downstream
MT2Expanded <- MT2_top
MT2Expanded$end <- MT2Expanded$start + 2e5 
MT2Expanded$start <- MT2Expanded$start - 2e5 
MT2ExpandedGR <- makeGRangesFromDataFrame(MT2Expanded, keep.extra.columns = T)

MT2ExpandedFilteredGRList <- as.list(MT2ExpandedGR)
names(MT2ExpandedFilteredGRList) <- MT2_top$repNames 

MT2ExpandedFilteredGRListStages <- rep(list(MT2ExpandedFilteredGRList), 5)
names(MT2ExpandedFilteredGRListStages) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")

# getting coverage counts
allMT2CoverageCounts_GV <- lapply(X = 1:length(MT2ExpandedFilteredGRListStages$"bam_GV"), 
                                  function(X) coverageForSummedPlot(X, coverage_GV, "bam_GV")) 
allMT2CoverageCounts_1C <- lapply(X = 1:length(MT2ExpandedFilteredGRListStages$"bam_1C"), 
                                  function(X) coverageForSummedPlot(X, coverage_1C, "bam_1C")) 
allMT2CoverageCounts_2C <- lapply(X = 1:length(MT2ExpandedFilteredGRListStages$"bam_2C"), 
                                  function(X) coverageForSummedPlot(X, coverage_2C, "bam_2C")) 
allMT2CoverageCounts_2C_aphi <- lapply(X = 1:length(MT2ExpandedFilteredGRListStages$"bam_2C_aphi"), 
                                       function(X) coverageForSummedPlot(X, coverage_2C_aphi, "bam_2C_aphi")) 
allMT2CoverageCounts_4C <- lapply(X = 1:length(MT2ExpandedFilteredGRListStages$"bam_4C"), 
                                  function(X) coverageForSummedPlot(X, coverage_4C, "bam_4C")) 

# all elements coverage by stages (long format)
all_stages_counts <- list(allMT2CoverageCounts_GV, allMT2CoverageCounts_1C, allMT2CoverageCounts_2C, allMT2CoverageCounts_2C_aphi, allMT2CoverageCounts_4C)
all_stages_names <- c("GV", "1C", "2C", "2C_aphi", "4C")
all_coverage_counts <- lapply(X = 1:length(all_stages_counts), 
                              function(X) positionCoverageDF(allMT2CoverageCounts = all_stages_counts[[X]], stage = all_stages_names[X]))

all_coverage_counts_df <- as.data.frame(do.call(rbind, all_coverage_counts))
rownames(all_coverage_counts_df) <- NULL

all_elements_coverage_list <- lapply(1:length(MT2OriginalGR), positionCoverageDFOneElementAllStages)
names(all_elements_coverage_list) <- paste0(MT2_top$seqnames, "\t",
                                            MT2_top$start + 1.5e5, "\t", 
                                            MT2_top$end - 1.5e5, "\t strand: ", 
                                            MT2_top$strand, "\t element: ", 
                                            MT2_top$name)

# plot
all_elements_plot_list <- lapply(1:length(MT2OriginalGR), singleElementPlotByStages)
pdf("single_elements_coverage_plot_no_filter1.pdf", width = 20, height = 10)
invisible(lapply(all_elements_plot_list, print))
dev.off()


