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
# options(bitmapType='cairo')

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/statistics")

# find coverage and plot function
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
  MT2SampleGR <- MT2_filtered_ranges_fpkm_filtered_grList_all[[bam_stage_name]][[x]]
  MT2SampleRL <- as(MT2SampleGR, "RangesList")
  MT2SampleCoverage <- Views(coverage_stage, MT2SampleRL)
  MT2SampleCoverage <- MT2SampleCoverage[[as.character(seqnames(MT2SampleGR[1]))]]
  MT2SampleCoverageCounts <- as.vector(unlist(viewApply(MT2SampleCoverage, as.vector)))
  MT2SampleCoveragePos <- unlist(apply(data.frame(x = start(MT2SampleCoverage), y = end(MT2SampleCoverage)), 1, function(x) x[1]:x[2]))
  
  # original sample data
  MT2SampleOriginal <- MT2_top[MT2_top$repNames == names(MT2_filtered_ranges_fpkm_filtered_grList_all[[bam_stage_name]][[x]])[1], ]
  MT2SampleExpanded <- MT2Expanded[MT2Expanded$repNames == names(MT2_filtered_ranges_fpkm_filtered_grList_all[[bam_stage_name]][[x]])[1], ]
  
  # counts with names = positions
  names(MT2SampleCoverageCounts) <- MT2SampleCoveragePos
  
  # data.frame for plot
  countsZero <- rep(0, 12e5 + 1)
  names(countsZero) <- MT2SampleExpanded$start : MT2SampleExpanded$end
  MT2SampleCoverageCountsFull <- replace(countsZero, names(countsZero) %in% names(MT2SampleCoverageCounts), MT2SampleCoverageCounts)
  
  MT2SampleCoveragePosition <- as.numeric(names(MT2SampleCoverageCountsFull))
  MT2SampleCoveragePosition <- MT2SampleCoveragePosition - MT2SampleExpanded$start - 6e5
  
  if (MT2SampleOriginal$strand == "-"){
    MT2SampleCoveragePosition <- MT2SampleCoveragePosition + (MT2SampleOriginal$end - MT2SampleOriginal$start)
    MT2SampleCoveragePosition <- rev(MT2SampleCoveragePosition)
  }
  
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
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-5e+05"):which(rownames(position_matrix) == as.character(500000 + element_width)), ]
  
  allMT2CoverageCountsSummedFinal <- rowSums(position_matrix)
  allMT2CoverageCountsSummedFinal <- data.frame(pos = as.numeric(names(allMT2CoverageCountsSummedFinal)), 
                                                count = allMT2CoverageCountsSummedFinal)
  allMT2CoverageCountsSummedFinal$element <- ifelse((allMT2CoverageCountsSummedFinal$pos >= 0 & allMT2CoverageCountsSummedFinal$pos <= element_width),
                                                    "in_element", 
                                                    "out_element")
  allMT2CoverageCountsSummedFinal$stage <- stage
  return(allMT2CoverageCountsSummedFinal)
}

plotCumFPKMBinArea <- function(tile_width){
  
  # upstream
  all_coverage_counts_df_merged_upstream <- all_coverage_counts_df_merged[all_coverage_counts_df_merged$pos >= -150000 & all_coverage_counts_df_merged$pos < 0, ]
  all_coverage_counts_df_merged_upstream <- 
    all_coverage_counts_df_merged_upstream %>%
    group_by(indx = gl(ceiling(nrow(all_coverage_counts_df_merged_upstream)/tile_width), tile_width, nrow(all_coverage_counts_df_merged_upstream))) %>%
    select(2:6) %>%
    summarise_each(funs(sum))
  all_coverage_counts_df_merged_upstream$position <- paste(seq(150000, tile_width, -tile_width), "upstream", sep = "_")
  
  # MT2
  all_coverage_counts_df_merged_mt2<- all_coverage_counts_df_merged[all_coverage_counts_df_merged$pos >= 0 & all_coverage_counts_df_merged$pos <= element_width, ]
  all_coverage_counts_df_merged_mt2 <- 
    all_coverage_counts_df_merged_mt2 %>%
    select(2:6) %>%
    summarise_each(funs(sum))
  all_coverage_counts_df_merged_mt2$position <- "element" 
  all_coverage_counts_df_merged_mt2$indx <- 1
  all_coverage_counts_df_merged_mt2 <- all_coverage_counts_df_merged_mt2[, c(7, 1:6)]
  
  # downstream
  all_coverage_counts_df_merged_downstream <- all_coverage_counts_df_merged[all_coverage_counts_df_merged$pos > element_width & all_coverage_counts_df_merged$pos <= (150000 + element_width), ]
  all_coverage_counts_df_merged_downstream <- 
    all_coverage_counts_df_merged_downstream %>%
    group_by(indx = gl(ceiling(nrow(all_coverage_counts_df_merged_downstream)/tile_width), tile_width, nrow(all_coverage_counts_df_merged_downstream))) %>%
    select(2:6) %>%
    summarise_each(funs(sum))
  all_coverage_counts_df_merged_downstream$position <- paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_")
  
  # creating one data.frame with mt2 and upstream/downstream positions
  all_coverage_counts_df_sum <- rbind(all_coverage_counts_df_merged_upstream, all_coverage_counts_df_merged_mt2, all_coverage_counts_df_merged_downstream)
  all_coverage_counts_df_sum <- all_coverage_counts_df_sum[, -1]
  
  # calculating FPKM for bins
  tile_kb <- tile_width / 1000
  all_coverage_counts_df_sum_fpkm <- all_coverage_counts_df_sum
  all_coverage_counts_df_sum_fpkm[, "GV"] <- all_coverage_counts_df_sum_fpkm[, "GV"] / (number_of_reads["GV"] * tile_kb)
  all_coverage_counts_df_sum_fpkm[, "1C"] <- all_coverage_counts_df_sum_fpkm[, "1C"] / (number_of_reads["1C"] * tile_kb)
  all_coverage_counts_df_sum_fpkm[, "2C"] <- all_coverage_counts_df_sum_fpkm[, "2C"]  / (number_of_reads["2C"] * tile_kb)
  all_coverage_counts_df_sum_fpkm[, "2C_aphi"] <- all_coverage_counts_df_sum_fpkm[, "2C_aphi"] / (number_of_reads["2C_aphi"] *  tile_kb)
  all_coverage_counts_df_sum_fpkm[, "4C"] <- all_coverage_counts_df_sum_fpkm[, "4C"] / (number_of_reads["4C"] * tile_kb)
  all_coverage_counts_df_sum_fpkm <- as.data.frame(all_coverage_counts_df_sum_fpkm)
  
  # melting data.frame to long format for ggplot
  all_coverage_counts_df_sum_fpkm_melt <- melt(all_coverage_counts_df_sum_fpkm, id.vars = "position")
  colnames(all_coverage_counts_df_sum_fpkm_melt) <- c("position", "stage", "fpkm")
  all_coverage_counts_df_sum_fpkm_melt$element <- ifelse((all_coverage_counts_df_sum_fpkm_melt$position == "element"), "in_element",  "out_element")
  all_coverage_counts_df_sum_fpkm_melt$width <- ifelse((all_coverage_counts_df_sum_fpkm_melt$position == "element"), element_width, tile_width)
  all_coverage_counts_df_sum_fpkm_melt$pos <- rep(c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width)), 5)
  all_coverage_counts_df_sum_fpkm_melt <- all_coverage_counts_df_sum_fpkm_melt[, c("pos", "fpkm", "element", "stage", "position", "width")]
  
  ## output .csv
  all_coverage_counts_df_sum_fpkm_output <- all_coverage_counts_df_sum_fpkm
  all_coverage_counts_df_sum_fpkm_output$pos <- c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width))
  write.csv(all_coverage_counts_df_sum_fpkm_output, paste0("./coverage/", element_name, "_cummulativeFPKM_bins", tile_width, ".csv"), row.names = F)
  
  # order by stage
  stage_order <- rep(c("2C_aphi", "2C", "1C", "4C", "GV"), each = (nrow(all_coverage_counts_df_sum_fpkm_melt) / 5))
  all_coverage_counts_df_sum_fpkm_melt$stage <- factor(all_coverage_counts_df_sum_fpkm_melt$stage, levels = stage_order) 
  all_coverage_counts_df_sum_fpkm_melt <- all_coverage_counts_df_sum_fpkm_melt[order(all_coverage_counts_df_sum_fpkm_melt$stage), ]  
  
  # # bin plot
  all_coverage_counts_df_sum_fpkm_melt$stage <- factor(all_coverage_counts_df_sum_fpkm_melt$stage, levels = c("GV", "4C", "1C", "2C", "2C_aphi"))
  ggplot(all_coverage_counts_df_sum_fpkm_melt, aes(x = pos, y = fpkm)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + width, ymax = fpkm, fill = stage)) +
    scale_fill_manual(values = c("black", "grey60", "grey30", "grey50", "orange"), 
                      breaks = c("GV", "1C", "2C", "2C_aphi", "4C")) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width), 
                       breaks = c(seq(-150000, -10000, 10000), 0, seq(element_width, (140000 + element_width), 10000)), 
                       labels = c(seq(150, 10, -10), "MuERV", seq(0, 140, 10)), 
                       name = "bin") + 
    scale_y_continuous(name = "Cummulative FPKM") +  
    coord_cartesian(ylim = c(0, 2000)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = paste0("./coverage/", element_name, "_cummulativeFPKM_bins", tile_width, "_binPlot.pdf"), width = 30, height = 10)
}

# # .bam files paths
# filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
#                          "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
#                          "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
#                          "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
#                          "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))
# 
# # # # reading .bam files
# bam_GV <- readGAlignmentPairs(filenames[1])
# bam_1C <- readGAlignmentPairs(filenames[2])
# bam_2C <- readGAlignmentPairs(filenames[3])
# bam_2C_aphi <- readGAlignmentPairs(filenames[4])
# bam_4C <- readGAlignmentPairs(filenames[5])
# # 
# # finding coverage of BAM files
# coverage_GV <- coverage(bam_GV)
# coverage_1C <- coverage(bam_1C)
# coverage_2C <- coverage(bam_2C)
# coverage_2C_aphi <- coverage(bam_2C_aphi)
# coverage_4C <- coverage(bam_4C)

################################################################################### reading full list, filtering by hand
# reading original list (top 100 MT2 and ORR1A0 solo LTRs by expression in 2Cell)
MT2_ORR1A0_all <- read.delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", stringsAsFactors = F, header = T)[, -4]
MT2_full <- read.delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_full_orderedByFPKMIn2cell.txt", stringsAsFactors = F, header = T)[, -c(4, 8)]
MT2_full$repName <- "MT2_full"

all_elements <- rbind(MT2_ORR1A0_all, MT2_full)
all_elements$fullName <- paste0(all_elements$seqnames, ":",
                                all_elements$start, "-", 
                                all_elements$end, ":", 
                                all_elements$strand, "|", 
                                all_elements$repName)

################################################################################### finding overlaps between original repeat ranges and expanded ranges
# expanding ranges by 150kb in both directions, making expanded and original GRanges
all_elements_expanded_gr <- all_elements
all_elements_expanded_gr$end <- all_elements_expanded_gr$start + 1.5e5 
all_elements_expanded_gr$start <- all_elements_expanded_gr$start - 1.5e5 
all_elements_expanded_gr <- makeGRangesFromDataFrame(all_elements_expanded_gr, keep.extra.columns = T)

all_elements_original_gr <- makeGRangesFromDataFrame(all_elements, keep.extra.columns = T)

# finding overlaps 
all_elements_overlaps <- findOverlaps(all_elements_expanded_gr, all_elements_original_gr)

# combining both in data.frame, removing regions overlaping with themself
all_elements_both <- as.data.frame(cbind(all_elements_original_gr[subjectHits(all_elements_overlaps)]$fullName, 
                                         all_elements_expanded_gr[queryHits(all_elements_overlaps)]$fullName))
colnames(all_elements_both) <- c("fullName_original", "fullName_expanded")
all_elements_both <- all_elements_both[all_elements_both$"fullName_original" != all_elements_both$"fullName_expanded", ]
all_elements_both_names <- as.character(unique(all_elements_both$fullName_expanded))

# filtering from original table all 
all_elements_filtered <- all_elements[!(all_elements$"fullName" %in% all_elements_both_names), ]

# names picked by hand for filtering
hand_picked_filter_elements_MT2_solo <- c("chr1:93965701-93966191:+|MT2_Mm", "chr2:157528079-157528349:+|MT2_Mm", 
                                          "chr5:137721893-137722008:-|MT2_Mm", "chr6:89202695-89202891:-|MT2_Mm", 
                                          "chr12:19108326-19108819:-|MT2_Mm", "chr12:19108326-19108819:-|MT2_Mm", 
                                          "chr1:83031132-83031678:+|MT2_Mm")
hand_picked_filter_elements_ORR1A0 <- c("chr17:6420849-6421203:+|ORR1A0", "chr17:6664850-6665198:-|ORR1A0", 
                                        "chr3:88577651-88577998:-|ORR1A0", "chr13:10040804-10041141:+|ORR1A0")
hand_picked_filter_elements_MT2_full <- c("chr11:60651277-60657678:-|MT2_full", "chr1:85172507-85178976:+|MT2_full", 
                                          "chr3:79046228-79052711:-|MT2_full", "chr12:19894744-19901175:+|MT2_full", 
                                          "chr13:76077210-76083662:-|MT2_full",  "chr13:114039987-114046421:+|MT2_full", 
                                          "chr13:119864904-119871300:+|MT2_full")

all_elements_filtered <- subset(all_elements_filtered, !(all_elements_filtered$fullName %in% c(hand_picked_filter_elements_MT2_solo, 
                                                                                               hand_picked_filter_elements_ORR1A0, 
                                                                                               hand_picked_filter_elements_MT2_full)))
# rm(list = ls()[grepl("hand_picked|both|overlaps|original|expanded|MT2|all_elements$", ls())])

# top 100 MT2/ORR1A0
all_elements_original <- 
  all_elements_filtered %>%
  group_by(repName) %>%
  arrange(desc(FPKM)) %>%
  slice(1:100) %>% 
  ungroup() %>% 
  select(-FPKM)

all_elements_original <- as.data.frame(all_elements_original)

# # making TxDb object from knownGene gtf from UCSC, getting FPKM
# knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz")
# knownGenes_gtf_gr <- genes(knownGenes_gtf)
# 
# # reading repeatMasker table, filtering for MT/MT2/ORR
# rptmsk <- read.delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", stringsAsFactors = F, header = T) 
# rptmsk <- rptmsk[, c("genoName", "genoStart", "genoEnd", "strand", "repName")]
# colnames(rptmsk) <- c("seqnames", "start", "end", "strand", "gene_id")
# rptmsk <- makeGRangesFromDataFrame(rptmsk, keep.extra.columns = T)
# 
# rptmskFiltered <- rptmsk[grep("MT|ORR", mcols(rptmsk)$gene_id)]
# rptmskFiltered <- rptmskFiltered[-grep("MMTV-int", mcols(rptmskFiltered)$gene_id)]
# 
# knownGenesAndRepeatMasker <- c(knownGenes_gtf_gr, rptmsk)
# 
##################################################################################### run from here
##################################################################################### 
##################################################################################### 
##################################################################################### 

# rm(list = ls()[!grepl("bam|coverage|count_function|filterRangesByFPKM|coverageForSummedPlot|positionCoverageDF|knownGenesAndRepeatMasker$|all_elements_original|plotCumFPKMBinArea|_MT2_full|_MT2_solo", ls())])
# rm(list = ls()[grepl("all_coverage_counts", ls())])

# reading original MT2 list (top 100 by expression in 2Cell)
element_name <- "MT2_full"

MT2_all <- all_elements_original[all_elements_original$repName == element_name, 1:5]
colnames(MT2_all) <- c("seqnames", "start", "end", "strand", "repNames")
MT2_all$repNames <- paste0(MT2_all$seqnames, ":",
                           MT2_all$start, "-", 
                           MT2_all$end, ":", 
                           MT2_all$strand)

# names picked by hand for filtering
MT2_hand_picked_filter <- c("chr11:60651277-60657678:-", "chr1:85172507-85178976:+", "chr3:79046228-79052711:-", 
                            "chr12:19894744-19901175:+", "chr12:19931732-19938183:-", "chr13:76077210-76083662:-")

MT2_filtered <- subset(MT2_all, !(MT2_all$repNames %in% MT2_hand_picked_filter))

# top 100
MT2_top <- MT2_filtered[1:100, ]
MT2OriginalGR <- makeGRangesFromDataFrame(MT2_top, keep.extra.columns = T)

# expand ranges 200kb up and downstream
MT2Expanded <- MT2_top
MT2Expanded$end <- MT2Expanded$start + 6e5 
MT2Expanded$start <- MT2Expanded$start - 6e5 
MT2ExpandedGR <- makeGRangesFromDataFrame(MT2Expanded, keep.extra.columns = T)

# joining filtered knownGenes and repeatMasker tables (to be used later for filtering reads)
knownGenesAndRepeatMaskerFiltered <- GenomicRanges::setdiff(knownGenesAndRepeatMasker, MT2OriginalGR, ignore.strand = T)
seqlevels(knownGenesAndRepeatMaskerFiltered, force = T) <- seqlevels(bam_2C)

# filtering with knownGenes and RepeatMasker
MT2ExpandedFilteredGRList <- lapply(MT2ExpandedGR, function(x) GenomicRanges::setdiff(x, knownGenesAndRepeatMaskerFiltered, ignore.strand = T))
names(MT2ExpandedFilteredGRList) <- MT2_top$repNames 
MT2ExpandedFilteredGRList <- MT2ExpandedFilteredGRList[!names(MT2ExpandedFilteredGRList) %in% MT2_hand_picked_filter]

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
element_width <- max(width(MT2OriginalGR))

all_stages_bam_names <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
MT2_filtered_ranges_fpkm_filtered_grList_all <- lapply(all_stages_bam_names, filterRangesByFPKM, fpkm_limit = 1)
names(MT2_filtered_ranges_fpkm_filtered_grList_all) <- all_stages_bam_names

# getting coverage counts
allMT2CoverageCounts_GV <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_GV"), 
                                  function(X) coverageForSummedPlot(X, coverage_GV, "bam_GV")) 
allMT2CoverageCounts_1C <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_1C"), 
                                  function(X) coverageForSummedPlot(X, coverage_1C, "bam_1C")) 
allMT2CoverageCounts_2C <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_2C"), 
                                  function(X) coverageForSummedPlot(X, coverage_2C, "bam_2C")) 
allMT2CoverageCounts_2C_aphi <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_2C_aphi"), 
                                       function(X) coverageForSummedPlot(X, coverage_2C_aphi, "bam_2C_aphi")) 
allMT2CoverageCounts_4C <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_4C"), 
                                  function(X) coverageForSummedPlot(X, coverage_4C, "bam_4C")) 

#########################################################################################################
# summing all stages coverage counts
all_stages_counts <- list(allMT2CoverageCounts_GV, allMT2CoverageCounts_1C, allMT2CoverageCounts_2C, allMT2CoverageCounts_2C_aphi, allMT2CoverageCounts_4C)
all_stages_names <- c("GV", "1C", "2C", "2C_aphi", "4C")
all_coverage_counts <- lapply(X = 1:length(all_stages_counts), 
                              function(X) positionCoverageDF(allMT2CoverageCounts = all_stages_counts[[X]], stage = all_stages_names[X]))
all_coverage_counts_df <- do.call(rbind, all_coverage_counts)
all_coverage_counts_df$stage <- factor(all_coverage_counts_df$stage, levels = c("GV", "1C", "2C", "2C_aphi", "4C"))
all_coverage_counts_df$count[all_coverage_counts_df$count == 0] <- NA
all_coverage_counts_df$pos[is.na(all_coverage_counts_df$count)] <- NA

# plotting with ggplot scaled
MT2SampleCoveragePlotSummed <- 
  ggplot(all_coverage_counts_df, aes(x = pos, y = count)) +
  geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
  scale_fill_manual(values = c("grey", "black")) +
  scale_x_continuous(limits = c(-5e5, 5e5 + element_width)) + 
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(stage ~ .)

ggsave(filename = paste0(element_name, "_500kExpanded_ylim100_summedCoverage.pdf"), 
       plot = MT2SampleCoveragePlotSummed, width = 30, height = 10)

##################################################################################################### bin plot
all_coverage_counts_list_merged <- lapply(X = 1:length(all_stages_counts), FUN = function(X) all_coverage_counts[[X]][, 1:2])
all_coverage_counts_df_merged <- Reduce(function(...) merge(..., by = 'pos', all.x = TRUE), all_coverage_counts_list_merged)
colnames(all_coverage_counts_df_merged) <- c("pos", "GV", "1C", "2C", "2C_aphi", "4C")
