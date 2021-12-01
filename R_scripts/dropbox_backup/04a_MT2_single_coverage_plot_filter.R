library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
library("ggplot2")
library("dplyr")
library("tidyr")
library("BiocParallel")
library("DESeq2")
library("reshape2")

rm(list = ls()[!grepl("bam|coverage", ls())])
rm(list = ls()[grepl("all|Summed|stage", ls())])

count_function <- function(bam_reads){
  register(MulticoreParam())
  se_filtered <- summarizeOverlaps(features = MT2_filtered_ranges_gr,
                                   reads = bam_reads, 
                                   mode = "Union", 
                                   singleEnd = FALSE, 
                                   ignore.strand = TRUE)
  return(assay(se_filtered))
}

filterRangesByFPKM <- function(bam_name, fpkm_limit){
  fpkm_df_limited_names <- MT2_filtered_ranges_df_fpkm_filtered[MT2_filtered_ranges_df_fpkm_filtered[, bam_name] > fpkm_limit, "feature_names"]
  MT2_filtered_ranges_fpkm_filtered_df <- MT2_filtered_ranges_df[!(rownames(MT2_filtered_ranges_df) %in% fpkm_df_limited_names), ]
  MT2_filtered_ranges_fpkm_filtered_gr <- makeGRangesFromDataFrame(MT2_filtered_ranges_fpkm_filtered_df)
  names(MT2_filtered_ranges_fpkm_filtered_gr) <- gsub("\\.[0-9]*", "", names(MT2_filtered_ranges_fpkm_filtered_gr))
  MT2_filtered_ranges_fpkm_filtered_grList <- split(MT2_filtered_ranges_fpkm_filtered_gr, names(MT2_filtered_ranges_fpkm_filtered_gr))
  return(MT2_filtered_ranges_fpkm_filtered_grList)
}

coverageForSummedPlot <- function(x, coverage_stage, bam_stage_name){
  MT2SampleGR <- MT2ExpandedFilteredGRListStages[[bam_stage_name]][[x]]
  seqlevels(MT2SampleGR, force = T) <- seqlevels(bam_2C)
  MT2SampleRL <- as(MT2SampleGR, "RangesList")
  MT2SampleCoverage <- Views(coverage_stage, MT2SampleRL)
  MT2SampleCoverage <- MT2SampleCoverage[[as.character(seqnames(MT2SampleGR[1]))]]
  MT2SampleCoverageCounts <- as.vector(unlist(viewApply(MT2SampleCoverage, as.vector)))
  MT2SampleCoveragePos <- unlist(apply(data.frame(x = start(MT2SampleCoverage), y = end(MT2SampleCoverage)), 1, function(x) x[1]:x[2]))
  
  # original sample data
  MT2SampleOriginal <- MT2_top[MT2_top$repNames == names(MT2ExpandedFilteredGRListStages[[bam_stage_name]][x])[1], ]
  MT2SampleExpanded <- MT2Expanded[MT2Expanded$repNames == names(MT2ExpandedFilteredGRListStages[[bam_stage_name]][x])[1], ]
  
  # counts with names = positions
  names(MT2SampleCoverageCounts) <- MT2SampleCoveragePos
  
  # data.frame for plot
  countsZero <- rep(0, 4e5 + 1)
  names(countsZero) <- MT2SampleExpanded$start : MT2SampleExpanded$end
  MT2SampleCoverageCountsFull <- replace(countsZero, names(countsZero) %in% names(MT2SampleCoverageCounts), MT2SampleCoverageCounts)
  
  MT2SampleCoveragePosition <- as.numeric(names(MT2SampleCoverageCountsFull))
  MT2SampleCoveragePosition <- MT2SampleCoveragePosition - MT2SampleExpanded$start - 2e5
  
  names(MT2SampleCoverageCountsFull) <- MT2SampleCoveragePosition
  return(MT2SampleCoverageCountsFull)
}

positionCoverageDF <- function(allMT2CoverageCounts, stage){
  # creating position matrix with all counts
  all_positions <- unique(unlist(lapply(allMT2CoverageCounts, names)))
  position_matrix <- matrix(NA, 
                            nrow = length(all_positions), 
                            ncol = length(allMT2CoverageCounts), 
                            dimnames = list(all_positions, c(1:length(allMT2CoverageCounts))))
  for (i in seq_along(allMT2CoverageCounts)) {
    position_matrix[names(allMT2CoverageCounts[[i]]), i] <- allMT2CoverageCounts[[i]]
  }
  position_matrix <- position_matrix[complete.cases(position_matrix), ]
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000"):which(rownames(position_matrix) == "150000"), ]
  return(position_matrix)
}

positionCoverageDFOneElementAllStages <- function(x){
  one_element <- all_coverage_counts_df[, x, drop = F]
  colnames(one_element) <- "count"
  one_element$pos <- rep(seq(-150000, 150000), 5)
  one_element$stage <- rep(all_stages_names, each = 300001)
  one_element$stage <- factor(one_element$stage, levels = c("GV", "1C", "2C", "2C_aphi", "4C"))
  one_element$element <- ifelse((one_element$pos >= 0 & one_element$pos <= width(MT2OriginalGR[x])), "in_element", "out_element")
  one_element$count[one_element$count == 0] <- NA
  one_element$pos[is.na(one_element$count)] <- NA
  return(one_element)
}

singleElementPlotByStages <- function(x){
  plot_df <- all_elements_coverage_list[[x]]
  single_element_facet_ggplot <- 
    ggplot(plot_df, aes(x = pos, y = count)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    scale_x_continuous(limits = c(-150000, 150000)) + 
    scale_fill_manual(values = c("red", "black")) +
    facet_grid(stage ~ .) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line.x = element_blank()) +
    coord_cartesian(ylim = c(0, 10)) +
    ggtitle(names(all_elements_coverage_list[x])) 
}

# reading ranges
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/coverage")
MT2_top <- read.csv("MT2_ORR_single_elements_range.csv", stringsAsFactors = F, header = T)
MT2OriginalGR <- makeGRangesFromDataFrame(MT2_top, keep.extra.columns = T)

# expand ranges 200kb up and downstream
MT2Expanded <- MT2_top
MT2Expanded$end <- MT2Expanded$start + 2e5 
MT2Expanded$start <- MT2Expanded$start - 2e5 
MT2ExpandedGR <- makeGRangesFromDataFrame(MT2Expanded, keep.extra.columns = T)

# making TxDb object from knownGene gtf from UCSC, getting FPKM
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mm10/UCSC_knownGenes_mm10_20160513.gtf.gz")
knownGenes_gtf_gr <- genes(knownGenes_gtf)

# reading repeatMasker table, filtering for MT/MT2/ORR
rptmsk <- import.bed("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20160607.bed.gz")
mcols(rptmsk) <- mcols(rptmsk)$name
names(mcols(rptmsk)) <- "gene_id"

# joining filtered knownGenes and repeatMasker tables (to be used later for filtering reads)
knownGenesAndRepeatMasker <- c(knownGenes_gtf_gr, rptmsk)
knownGenesAndRepeatMaskerFiltered <- GenomicRanges::setdiff(knownGenesAndRepeatMasker, MT2OriginalGR, ignore.strand = T)
seqlevels(knownGenesAndRepeatMaskerFiltered, force = T) <- seqlevels(bam_2C)

# filtering with knownGenes and RepeatMasker
MT2ExpandedFilteredGRList <- lapply(MT2ExpandedGR, function(x) GenomicRanges::setdiff(x, knownGenesAndRepeatMaskerFiltered, ignore.strand = T))
names(MT2ExpandedFilteredGRList) <- MT2_top$repNames 

# getting FPKM values for filtered features
MT2_filtered_ranges <- lapply(MT2ExpandedFilteredGRList, as.data.frame)
MT2_filtered_ranges_df <- do.call(rbind, MT2_filtered_ranges)
MT2_filtered_ranges_gr <- makeGRangesFromDataFrame(MT2_filtered_ranges_df)

### calculating fpkm values
# getting library size in millions of reads
logs_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out", 
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAmLog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WELog.final.out"))

number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- c("GV", "1C", "2C", "2C_aphi", "4C")
number_of_reads <- number_of_reads / 10^6

# getting counts of all samples with summarizeOverlaps function
count_list <- lapply(list(bam_GV, bam_1C, bam_2C, bam_2C_aphi, bam_4C), count_function)
count_df <- do.call(cbind, count_list)
colnames(count_df) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
count_df <- as.data.frame(count_df)
count_df$feature_names <- rownames(count_df)
rownames(count_df) <- NULL
MT2_filtered_ranges_df_counts <- cbind(MT2_filtered_ranges_df, count_df)
rownames(MT2_filtered_ranges_df_counts) <- NULL

# calculating fpkm from counts
MT2_filtered_ranges_df_fpkm <- MT2_filtered_ranges_df_counts
MT2_filtered_ranges_df_fpkm[, "bam_GV"] <- MT2_filtered_ranges_df_fpkm$bam_GV / 
  (number_of_reads["GV"] * (MT2_filtered_ranges_df_fpkm$width / 1000))
MT2_filtered_ranges_df_fpkm[, "bam_1C"] <- MT2_filtered_ranges_df_fpkm$bam_1C / 
  (number_of_reads["1C"] * (MT2_filtered_ranges_df_fpkm$width / 1000))
MT2_filtered_ranges_df_fpkm[, "bam_2C"] <- MT2_filtered_ranges_df_fpkm$bam_2C / 
  (number_of_reads["2C"] * (MT2_filtered_ranges_df_fpkm$width / 1000))
MT2_filtered_ranges_df_fpkm[, "bam_2C_aphi"] <- MT2_filtered_ranges_df_fpkm$bam_2C_aphi / 
  (number_of_reads["2C_aphi"] * (MT2_filtered_ranges_df_fpkm$width / 1000))
MT2_filtered_ranges_df_fpkm[, "bam_4C"] <- MT2_filtered_ranges_df_fpkm$bam_4C / 
  (number_of_reads["4C"] * (MT2_filtered_ranges_df_fpkm$width / 1000))

# removing MT2 positions from fpkm 
MT2_filtered_ranges_df_fpkm_gr <- makeGRangesFromDataFrame(MT2_filtered_ranges_df_fpkm, keep.extra.columns = T)
MT2_filtered_ranges_df_fpkm_gr_overlaps <- findOverlaps(MT2_filtered_ranges_df_fpkm_gr, MT2OriginalGR)
MT2_filtered_ranges_df_fpkm_gr <- MT2_filtered_ranges_df_fpkm_gr[-queryHits(MT2_filtered_ranges_df_fpkm_gr_overlaps)]
MT2_filtered_ranges_df_fpkm_filtered <- MT2_filtered_ranges_df_fpkm[MT2_filtered_ranges_df_fpkm$feature_names %in% MT2_filtered_ranges_df_fpkm_gr$feature_names, ]

#########################################################################################################
# filtering ranges by FPKM values
all_stages_bam_names <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
MT2_filtered_ranges_fpkm_filtered_grList_all <- lapply(all_stages_bam_names, filterRangesByFPKM, fpkm_limit = 1)
names(MT2_filtered_ranges_fpkm_filtered_grList_all) <- all_stages_bam_names
MT2ExpandedFilteredGRListStages <- MT2_filtered_ranges_fpkm_filtered_grList_all

# # no FPKM filtration
# MT2ExpandedFilteredGRListStages <- rep(list(MT2ExpandedFilteredGRList), 5)
# names(MT2ExpandedFilteredGRListStages) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")

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

# list of elements with data.frames by stages (long format)
all_elements_coverage_list <- lapply(1:length(MT2OriginalGR), positionCoverageDFOneElementAllStages)
names(all_elements_coverage_list) <- paste0(MT2_top$seqnames, "\t",
                                            MT2_top$start + 1.5e5, "\t", 
                                            MT2_top$end - 1.5e5, "\t strand: ", 
                                            MT2_top$strand, "\t element: ", 
                                            MT2_top$name)
# plot
all_elements_plot_list <- lapply(1:length(MT2OriginalGR), singleElementPlotByStages)
pdf("single_elements_coverage_plot_filter_fpkm1.pdf", width = 20, height = 10)
invisible(lapply(all_elements_plot_list, print))
dev.off()
