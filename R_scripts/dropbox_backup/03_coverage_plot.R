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

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/coverage")
rm(list = ls()[!grepl("bam|coverage", ls())])
rm(list = ls()[grepl("all|Summed|stage|bam_filenames|bamfiles|MT2_ORR1A0|bam_name", ls())])

# find coverage and plot function
countFunction <- function(bam_reads, features_gr){
  register(MulticoreParam())
  se_filtered <- summarizeOverlaps(features = features_gr,
                                   reads = bam_reads, 
                                   mode = "Union", 
                                   singleEnd = FALSE, 
                                   ignore.strand = TRUE)
  return(assay(se_filtered))
}

filterRangesByFPKM <- function(bam_name, fpkm_df, fpkm_limit, original_df){
  fpkm_df_limited_names <- fpkm_df[fpkm_df[, bam_name] > fpkm_limit, "feature_names"]
  fpkm_df_filtered <- original_df[!(rownames(original_df) %in% fpkm_df_limited_names), ]
  fpkm_df_filtered_gr <- makeGRangesFromDataFrame(fpkm_df_filtered)
  names(fpkm_df_filtered_gr) <- gsub("\\.[0-9]*", "", names(fpkm_df_filtered_gr))
  fpkm_df_filtered_grList <- split(fpkm_df_filtered_gr, names(fpkm_df_filtered_gr))
  return(fpkm_df_filtered_grList)
}

coverageForSummedPlot <- function(x, coverage_stage, bam_stage_name, granges_list){
  sample_gr <- granges_list[[bam_stage_name]][[x]]
  sample_rl <- as(sample_gr, "RangesList")
  sample_coverage <- Views(coverage_stage, sample_rl)
  sample_coverage <- sample_coverage[[as.character(seqnames(sample_gr[1]))]]
  sample_coverage_counts <- as.vector(unlist(viewApply(sample_coverage, as.vector)))
  sample_coverage_position <- unlist(apply(data.frame(x = start(sample_coverage), y = end(sample_coverage)), 1, function(x) x[1]:x[2]))
  
  # original sample data
  sample_range_original <- all_elements_original[all_elements_original$fullName == names(granges_list[[bam_stage_name]][[x]])[1], ]
  sample_range_expanded <- all_elements_expanded[all_elements_original$fullName == names(granges_list[[bam_stage_name]][[x]])[1], ]
  
  # counts with names = positions
  names(sample_coverage_counts) <- sample_coverage_position
  
  # data.frame for plot
  counts_zero <- rep(0, 4e5 + 1)
  names(counts_zero) <- sample_range_expanded$start : sample_range_expanded$end
  sample_coverage_counts_full <- replace(counts_zero, names(counts_zero) %in% names(sample_coverage_counts), sample_coverage_counts)
  
  sample_coverage_position <- as.numeric(names(sample_coverage_counts_full))
  sample_coverage_position <- sample_coverage_position - sample_range_expanded$start - 2e5
  
  names(sample_coverage_counts_full) <- sample_coverage_position
  return(sample_coverage_counts_full)
}

positionCoverageDF <- function(coverage_counts, stage, element_width){
  # creating position matrix with all counts
  all_positions <- unique(unlist(lapply(coverage_counts, names)))
  position_matrix <- matrix(NA, 
                            nrow = length(all_positions), 
                            ncol = length(coverage_counts), 
                            dimnames = list(all_positions, c(1:length(coverage_counts))))
  for (i in seq_along(coverage_counts)) {
    position_matrix[names(coverage_counts[[i]]), i] <- coverage_counts[[i]]
  }
  position_matrix <- position_matrix[complete.cases(position_matrix), ]
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000"):which(rownames(position_matrix) == as.character(150000 + element_width)), ]
  
  coverage_counts_summed <- rowSums(position_matrix)
  coverage_counts_summed <- data.frame(pos = as.numeric(names(coverage_counts_summed)), count = coverage_counts_summed)
  coverage_counts_summed$element <- ifelse((coverage_counts_summed$pos >= 0 & coverage_counts_summed$pos <= element_width), "in_element", "out_element")
  coverage_counts_summed$stage <- stage
  return(coverage_counts_summed)
}

# .bam files paths
filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))

# reading .bam files
bam_GV <- readGAlignmentPairs(filenames[1])
bam_1C <- readGAlignmentPairs(filenames[2])
bam_2C <- readGAlignmentPairs(filenames[3])
bam_2C_aphi <- readGAlignmentPairs(filenames[4])
bam_4C <- readGAlignmentPairs(filenames[5])

# finding coverage of BAM files
coverage_GV <- coverage(bam_GV)
coverage_1C <- coverage(bam_1C)
coverage_2C <- coverage(bam_2C)
coverage_2C_aphi <- coverage(bam_2C_aphi)
coverage_4C <- coverage(bam_4C)

################################################################################### reading full list, filtering by hand
# reading original list (top 100 MT2 and ORR1A0 solo LTRs by expression in 2Cell)
MT2_ORR1A0_all <- read.delim("MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", stringsAsFactors = F, header = T)[, -4]
MT2_full <- read.delim("MT2_full_orderedByFPKMIn2cell.txt", stringsAsFactors = F, header = T)[, -c(4, 8)]
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

################################################################################### filering by name
# names picked by hand for filtering
hand_picked_filter_elements_MT2_solo <- c("chr1:93965701-93966191:+|MT2_Mm", "chr2:157528079-157528349:+|MT2_Mm", 
                                          "chr5:137721893-137722008:-|MT2_Mm", "chr6:89202695-89202891:-|MT2_Mm", 
                                          "chr12:19108326-19108819:-|MT2_Mm", "chr12:19108326-19108819:-|MT2_Mm", 
                                          "chr1:83031132-83031678:+|MT2_Mm")
hand_picked_filter_elements_ORR1A0 <- c("chr17:6420849-6421203:+|ORR1A0", "chr17:6664850-6665198:-|ORR1A0", 
                                        "chr3:88577651-88577998:-|ORR1A0", "chr13:10040804-10041141:+|ORR1A0")
hand_picked_filter_elements_MT2_full <- c("chr11:60651277-60657678:-|MT2_full", "chr1:85172507-85178976:+|MT2_full", 
                                          "chr3:79046228-79052711:-|MT2_full", "chr12:19894744-19901175:+|MT2_full", 
                                          "chr13:76077210-76083662:-|MT2_full")
  
all_elements_filtered <- subset(all_elements_filtered, !(all_elements_filtered$fullName %in% c(hand_picked_filter_elements_MT2_solo, 
                                                                                               hand_picked_filter_elements_ORR1A0, 
                                                                                               hand_picked_filter_elements_MT2_full)))
rm(list = ls()[grepl("hand_picked|both|overlaps|original|expanded|MT2|all_elements$", ls())])

################################################################################### creating original list of 100 MT2 and 100 ORR1A0 elements
# top 100 MT2/ORR1A0
all_elements_original <- 
  all_elements_filtered %>%
  group_by(repName) %>%
  arrange(desc(FPKM)) %>%
  slice(1:100) %>% 
  ungroup() %>% 
  select(-FPKM)

all_elements_original <- as.data.frame(all_elements_original)
all_elements_original_gr <- makeGRangesFromDataFrame(all_elements_original, keep.extra.columns = T)

# expand ranges 200kb up and downstream
all_elements_expanded <- all_elements_original
all_elements_expanded$end <- all_elements_expanded$start + 2e5 
all_elements_expanded$start <- all_elements_expanded$start - 2e5 
all_elements_expanded_gr <- makeGRangesFromDataFrame(all_elements_expanded, keep.extra.columns = T)

################################################################################### filtering by knownGene and repeatMasker tables
# making TxDb object from knownGene gtf from UCSC, getting FPKM
knownGenes_gtf_gr <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mm10/UCSC_knownGenes_mm10_20160513.gtf.gz")
knownGenes_gtf_gr <- genes(knownGenes_gtf_gr)

# reading repeatMasker table
rptmsk_gr <- import.bed("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20160607.bed.gz")
mcols(rptmsk_gr) <- mcols(rptmsk_gr)$name
names(mcols(rptmsk_gr)) <- "gene_id"

# joining filtered knownGenes and repeatMasker tables (to be used later for filtering reads)
knownGenes_RepeatMasker_gr <- c(knownGenes_gtf_gr, rptmsk_gr)
knownGenes_RepeatMasker_gr <- GenomicRanges::setdiff(knownGenes_RepeatMasker_gr, all_elements_original_gr, ignore.strand = T)
seqlevels(knownGenes_RepeatMasker_gr, force = T) <- seqlevels(bam_2C)

# filtering with knownGenes and RepeatMasker
all_elements_expanded_filtered <- lapply(all_elements_expanded_gr, function(x) GenomicRanges::setdiff(x, knownGenes_RepeatMasker_gr, ignore.strand = T))
names(all_elements_expanded_filtered) <- all_elements_original$fullName 

################################################################################### filtering by FPKM values
# making GRanges from GRangesList
all_elements_expanded_filtered_gr <- lapply(all_elements_expanded_filtered, as.data.frame)
all_elements_expanded_filtered_df <- do.call(rbind, all_elements_expanded_filtered_gr)
all_elements_expanded_filtered_gr <- makeGRangesFromDataFrame(all_elements_expanded_filtered_df)

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
count_list <- lapply(list(bam_GV, bam_1C, bam_2C, bam_2C_aphi, bam_4C), countFunction, features_gr = all_elements_expanded_filtered_gr)
count_df <- do.call(cbind, count_list)
colnames(count_df) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
count_df <- as.data.frame(count_df)
count_df$feature_names <- rownames(count_df)
rownames(count_df) <- NULL

# binding counts with MT2/ORR1A0 ranges
all_elements_fpkm <- cbind(as.data.frame(all_elements_expanded_filtered_gr), count_df)
rownames(all_elements_fpkm) <- NULL

# calculating fpkm from counts
all_elements_fpkm[, "bam_GV"] <- all_elements_fpkm$bam_GV / (number_of_reads["GV"] * (all_elements_fpkm$width / 1000))
all_elements_fpkm[, "bam_1C"] <- all_elements_fpkm$bam_1C / (number_of_reads["1C"] * (all_elements_fpkm$width / 1000))
all_elements_fpkm[, "bam_2C"] <- all_elements_fpkm$bam_2C / (number_of_reads["2C"] * (all_elements_fpkm$width / 1000))
all_elements_fpkm[, "bam_2C_aphi"] <- all_elements_fpkm$bam_2C_aphi / (number_of_reads["2C_aphi"] * (all_elements_fpkm$width / 1000))
all_elements_fpkm[, "bam_4C"] <- all_elements_fpkm$bam_4C / (number_of_reads["4C"] * (all_elements_fpkm$width / 1000))

# removing repeat elements positions from FPKM 
all_elements_fpkm_gr <- makeGRangesFromDataFrame(all_elements_fpkm, keep.extra.columns = T)
all_elements_fpkm_gr_overlaps <- findOverlaps(all_elements_fpkm_gr, all_elements_original_gr)
all_elements_fpkm_gr <- all_elements_fpkm_gr[-queryHits(all_elements_fpkm_gr_overlaps)]
all_elements_fpkm_filtered <- as.data.frame(all_elements_fpkm_gr)

# filtering ranges by FPKM values
all_stages_bam_names <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
all_elements_GRList_FPKM_filtered <- lapply(all_stages_bam_names, 
                                            FUN = filterRangesByFPKM, 
                                            fpkm_limit = 1, 
                                            fpkm_df = all_elements_fpkm_filtered, 
                                            original_df = all_elements_expanded_filtered_df)
names(all_elements_GRList_FPKM_filtered) <- all_stages_bam_names

rm(list = ls()[grepl("knownGenes|rptmsk|all_elements_expanded_filtered|logs|number_of_reads|count_|all_elements_fpkm|all_stages_bam_names|all_elements_filtered", ls())])

#########################################################################################################
# getting coverage counts
all_coverage_GV <- lapply(X = 1:length(all_elements_GRList_FPKM_filtered$"bam_GV"), 
                          function(X) coverageForSummedPlot(X, coverage_stage = coverage_GV,
                                                            bam_stage_name = "bam_GV",
                                                            granges_list = all_elements_GRList_FPKM_filtered)) 
  
all_coverage_1C <- lapply(X = 1:length(all_elements_GRList_FPKM_filtered$"bam_1C"), 
                          function(X) coverageForSummedPlot(X, coverage_stage = coverage_1C,
                                                            bam_stage_name = "bam_1C",
                                                            granges_list = all_elements_GRList_FPKM_filtered)) 

all_coverage_2C <- lapply(X = 1:length(all_elements_GRList_FPKM_filtered$"bam_2C"), 
                          function(X) coverageForSummedPlot(X, coverage_stage = coverage_2C, 
                                                            bam_stage_name = "bam_2C", 
                                                            granges_list = all_elements_GRList_FPKM_filtered)) 

all_coverage_2C_aphi <- lapply(X = 1:length(all_elements_GRList_FPKM_filtered$"bam_2C_aphi"), 
                               function(X) coverageForSummedPlot(X, coverage_stage = coverage_2C_aphi, 
                                                                 bam_stage_name = "bam_2C_aphi", 
                                                                 granges_list = all_elements_GRList_FPKM_filtered)) 

all_coverage_4C <- lapply(X = 1:length(all_elements_GRList_FPKM_filtered$"bam_4C"), 
                          function(X) coverageForSummedPlot(X, coverage_stage = coverage_4C, 
                                                            bam_stage_name = "bam_4C", 
                                                            granges_list = all_elements_GRList_FPKM_filtered)) 

names(all_coverage_GV) <- names(all_elements_GRList_FPKM_filtered$"bam_GV")
names(all_coverage_1C) <- names(all_elements_GRList_FPKM_filtered$"bam_1C")
names(all_coverage_2C) <- names(all_elements_GRList_FPKM_filtered$"bam_2C")
names(all_coverage_2C_aphi) <- names(all_elements_GRList_FPKM_filtered$"bam_2C_aphi")
names(all_coverage_4C) <- names(all_elements_GRList_FPKM_filtered$"bam_4C")

#########################################################################################################
rm(list = ls()[grepl("all_elements_GRList_FPKM_filtered", ls())])

# summing all stages coverage counts
positionCoverageDFSum <- function(repeat_name){
  all_stages_counts <- list(all_coverage_GV[grep(repeat_name, names(all_coverage_GV))], 
                            all_coverage_1C[grep(repeat_name, names(all_coverage_1C))], 
                            all_coverage_2C[grep(repeat_name, names(all_coverage_2C))], 
                            all_coverage_2C_aphi[grep(repeat_name, names(all_coverage_2C_aphi))], 
                            all_coverage_4C[grep(repeat_name, names(all_coverage_4C))])
  all_stages_names <- c("GV", "1C", "2C", "2C_aphi", "4C")
  all_coverage_counts <- lapply(X = 1:length(all_stages_counts), 
                                function(X) positionCoverageDF(coverage_counts = all_stages_counts[[X]], 
                                                               stage = all_stages_names[X], 
                                                               element_width = max(width(all_elements_original_gr[grep(repeat_name, all_elements_original_gr$fullName)]))))
  return(all_coverage_counts)
}

MT2_solo_all_coverage_counts <- positionCoverageDFSum("MT2_Mm")
ORR1A0_all_coverage_counts <- positionCoverageDFSum("ORR1A0")
MT2_full_all_coverage_counts <- positionCoverageDFSum("MT2_full")

######################################################################################################### making data.frames for each element and plot coverage
dataFramePlotCoverage <- function(coverage_counts, repeat_name){
  
  # making data.frame for ggplot
  all_coverage_counts_df <- do.call(rbind, coverage_counts)
  all_coverage_counts_df$stage <- factor(all_coverage_counts_df$stage, levels = c("GV", "1C", "2C", "2C_aphi", "4C"))
  all_coverage_counts_df$count[all_coverage_counts_df$count == 0] <- NA
  all_coverage_counts_df$pos[is.na(all_coverage_counts_df$count)] <- NA
  
  # plotting summed coverage faceted
  ggplot(all_coverage_counts_df, aes(x = pos, y = count)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    scale_fill_manual(values = c("grey", "black")) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + max(width(all_elements_original_gr[grep(repeat_name, all_elements_original_gr$fullName)])))) + 
    coord_cartesian(ylim = c(0, 50)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_grid(stage ~ .) + 
    ggsave(filename = paste0(repeat_name, "_top100_ylim_50_fpkm_1_filtered_summed_coverage.pdf"), width = 30, height = 10)
  return(NULL)
}

dataFramePlotCoverage(MT2_solo_all_coverage_counts, "MT2_Mm")
dataFramePlotCoverage(ORR1A0_all_coverage_counts, "ORR1A0")
dataFramePlotCoverage(MT2_full_all_coverage_counts, "MT2_full")

rm(list = ls()[grepl("_all_coverage_counts", ls())])
