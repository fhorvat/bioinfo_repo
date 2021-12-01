library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
library("ggplot2")
library("dplyr")
library("tidyr")
library("reshape2")
library("BiocParallel")

options(bitmapType = "cairo")
rm(list = ls()[!grepl("^bam_|^coverage_|^repeatMasker$|^knownGenes|countFunction|filterRangesByFPKM|coverageForSummedPlot|positionCoverageDF|plotCummulativeFPKMBin", ls())])

rm(list = ls()[!grepl("^bam_|^coverage_", ls())])

################################################################################### functions
countFunction <- function(feature, bam_reads){
  register(MulticoreParam())
  se_filtered <- summarizeOverlaps(features = feature,
                                   reads = bam_reads, 
                                   mode = "Union", 
                                   singleEnd = FALSE, 
                                   ignore.strand = TRUE)
  return(assay(se_filtered))
}

filterRangesByFPKM <- function(bam_name, element_fpkm, element_ranges, fpkm_limit){
  fpkm_df_names_limit <- element_fpkm[element_fpkm[, bam_name] > fpkm_limit, "feature_names"]
  element_ranges_fpkm_filtered <- element_ranges[!(rownames(element_ranges) %in% fpkm_df_names_limit), ]
  element_ranges_fpkm_filtered_gr <- makeGRangesFromDataFrame(element_ranges_fpkm_filtered)
  names(element_ranges_fpkm_filtered_gr) <- gsub("\\.[0-9]*", "", names(element_ranges_fpkm_filtered_gr))
  element_ranges_fpkm_filtered_grList <- split(element_ranges_fpkm_filtered_gr, names(element_ranges_fpkm_filtered_gr))
  return(element_ranges_fpkm_filtered_grList)
}

coverageForSummedPlot <- function(x, coverage_stage, bam_stage_name, element_grList){
  
  # function takes (filtered) ranges of one element and returns relative 
  # coverage on each position of expanded range (+-2e5 kb) in .bam file
  
  # take element ranges, get original element data
  element_gr <- element_grList[[bam_stage_name]][[x]]
  element_original <- element_original_ranges[element_original_ranges$fullName == names(element_gr)[1]]
  element_expanded <- element_expanded_ranges[element_original_ranges$fullName == names(element_gr)[1]]
  
  # find coverage
  element_coverage <- Views(coverage_stage, as(element_gr, "RangesList"))
  element_coverage <- element_coverage[[as.character(seqnames(element_gr[1]))]]
  
  # separate coverage counts and positions
  element_coverage_counts <- as.vector(unlist(viewApply(element_coverage, as.vector)))
  element_coverage_pos <- unlist(apply(data.frame(x = start(element_coverage), y = end(element_coverage)), 1, function(x) x[1]:x[2]))
  names(element_coverage_counts) <- element_coverage_pos
  
  ### make coverage positions with width = +- 2e5 kb with 0 where there is no coverage
  counts_zero <- rep(0, 4e5 + 1)
  names(counts_zero) <- start(element_expanded):end(element_expanded)
  
  # replace 0 counts with real counts based on the name (= position in expanded element), get coverage position from names
  element_coverage_counts_full <- replace(counts_zero, names(counts_zero) %in% names(element_coverage_counts), element_coverage_counts)
  element_coverage_position_full <- as.numeric(names(element_coverage_counts_full))
  
  # make coverage position relative with element start = 0 and 2e5 kb up-/down-stream  
  element_coverage_position_full <- element_coverage_position_full - start(element_expanded) - 2e5
  
  # if element is on minus strand shift position for width of element and reverse coverage position 
  if (as.character(strand(element_original)) == "-"){
    element_coverage_position_full <- element_coverage_position_full + (end(element_original) - start(element_original))
    element_coverage_position_full <- rev(element_coverage_position_full)
  }
  
  # return vector with coverage counts and position as names
  names(element_coverage_counts_full) <- element_coverage_position_full
  return(element_coverage_counts_full)
}

positionCoverageDF <- function(element_coverage){
  
  # getting all unique position across all elements in one stage
  all_positions <- unique(unlist(lapply(element_coverage, names)))
  
  # creating position matrix for all counts (rows = position, columns = elements)
  position_matrix <- matrix(NA, 
                            nrow = length(all_positions), 
                            ncol = length(element_coverage), 
                            dimnames = list(all_positions, c(1:length(element_coverage))))
  
  # filling position matrix with coverage based on position names
  for (i in seq_along(element_coverage)) {
    position_matrix[names(element_coverage[[i]]), i] <- element_coverage[[i]]
  }
  
  # keep only position for which all elements have position and then keep only +-1.5e5 kb
  position_matrix <- position_matrix[complete.cases(position_matrix), ]
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000"):which(rownames(position_matrix) == as.character(150000 + element_width)), ]
  
  # sum rows = summed coverage of all elements per position
  element_coverage_summed <- rowSums(position_matrix)
  element_coverage_summed <- data.frame(pos = as.numeric(names(element_coverage_summed)), count = element_coverage_summed)
  element_coverage_summed$element <- ifelse(element_coverage_summed$pos >= 0 & element_coverage_summed$pos <= element_width, "in_element", "out_element")
  
  return(element_coverage_summed)
}

plotCummulativeFPKMBin <- function(tile_width, merged_coverage_df){
  
  all_coverage_counts_df_sum <- 
    rbind(merged_coverage_df %>% #upstream
            filter(pos >= -150000 & pos < 0) %>%
            group_by(indx = gl(ceiling(n()/tile_width), tile_width, n())) %>%
            select(2:6) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = paste(seq(150000, tile_width, -tile_width), "upstream", sep = "_")),
          merged_coverage_df %>% #element
            filter(pos >= 0 & pos <= element_width) %>%
            select(2:6) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = "element", 
                   indx = 1) %>%
            select(c(7, 1:6)), 
          merged_coverage_df %>% #downstream
            filter(pos > element_width & pos <= (150000 + element_width)) %>%
            group_by(indx = gl(ceiling(n()/tile_width), tile_width, n())) %>%
            select(2:6) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_"))) %>%
    select(-1)
  
  # calculating FPKM for bins
  invisible(lapply(X = names(number_of_reads),
                   FUN = function(X) all_coverage_counts_df_sum[, X] <<- all_coverage_counts_df_sum[, X] / (number_of_reads[X] * (tile_width / 1000))))
  
  # output .csv
  all_coverage_counts_df_sum %>%
    mutate(pos = c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width))) %>%
    write.csv(file = paste0(file_name, "_cummulativeFPKM_bins", tile_width, ".csv"), row.names = F)
  
  # melting data.frame to long format for ggplot
  all_coverage_counts_df_sum <- 
    melt(all_coverage_counts_df_sum, id.vars = "position") %>%
    select(position, stage = variable, fpkm = value) %>%
    mutate(element = ifelse((position == "element"), "in_element",  "out_element"), 
           width = ifelse((position == "element"), element_width, tile_width), 
           pos = rep(c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width)), 5), 
           stage = gsub("bam_", "", stage)) %>%
    select(pos, fpkm, element, stage, position, width) %>%
    mutate(stage = factor(stage, levels = c("2C_aphi", "2C", "1C", "4C", "GV"))) %>%
    arrange(stage) %>%
    mutate(stage = factor(stage, levels = c("GV", "4C", "1C", "2C", "2C_aphi")))
  
  # bin plot
  ggplot(all_coverage_counts_df_sum, aes(x = pos, y = fpkm)) +
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
    ggsave(filename = paste0(file_name, "_cummulativeFPKM_bins", tile_width, ".pdf"), width = 30, height = 10)
}

################################################################################## reading data
# .bam files paths
filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))

# .bam files
bam_GV <- readGAlignmentPairs(filenames[1])
bam_1C <- readGAlignmentPairs(filenames[2])
bam_2C <- readGAlignmentPairs(filenames[3])
bam_2C_aphi <- readGAlignmentPairs(filenames[4])
bam_4C <- readGAlignmentPairs(filenames[5])

# coverage of BAM files
coverage_GV <- coverage(bam_GV)
coverage_1C <- coverage(bam_1C)
coverage_2C <- coverage(bam_2C)
coverage_2C_aphi <- coverage(bam_2C_aphi)
coverage_4C <- coverage(bam_4C)

# repeatMasker table
repeatMasker <- read.delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", header = T, stringsAsFactors = F)
repeatMasker <- makeGRangesFromDataFrame(repeatMasker, keep.extra.columns = T)
mcols(repeatMasker) <- mcols(repeatMasker)$element_name
names(mcols(repeatMasker)) <- "element_name"

# making TxDb object from knownGene gtf from UCSC, getting exons and genes
knownGenesTxDb <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz")
knownGenesExons <- unlist(exonsBy(knownGenesTxDb, by = "gene"))
mcols(knownGenesExons) <- names(knownGenesExons)
names(mcols(knownGenesExons)) <- "element_name"
knownGenes <- genes(knownGenesTxDb)
names(mcols(knownGenes)) <- "element_name"

# getting library size in millions of reads
logs_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out", 
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAmLog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WELog.final.out"))

number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
number_of_reads <- number_of_reads / 10^6

# reading original element lists (MT2 and ORR1A0 solo LTRs ordered by expression in 2Cell)
MT2_ORR1A0_all <- read.delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", stringsAsFactors = F, header = T)[, -4]
MT2_full <- read.delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_full_orderedByFPKMIn2cell.txt", stringsAsFactors = F, header = T)[, -c(4, 8)]
MT2_full$repName <- "MT2_full"

################################################################################### all elements ranges and expanded ranges
# table with all elements (MT2/ORR1A solo LTRs, full MT2 elements)
all_elements <- rbind(MT2_ORR1A0_all, MT2_full)
all_elements$fullName <- paste0(all_elements$seqnames, ":",
                                all_elements$start, "-", 
                                all_elements$end, ":", 
                                all_elements$strand, "|", 
                                all_elements$repName)
all_elements <- all_elements[!grepl("random|Un", all_elements$seqnames), ]

# expanding ranges by 150kb in both directions, making expanded GRanges
all_elements_expanded_gr <- all_elements
all_elements_expanded_gr$end <- all_elements_expanded_gr$start + 1.5e5 
all_elements_expanded_gr$start <- all_elements_expanded_gr$start - 1.5e5 
all_elements_expanded_gr <- makeGRangesFromDataFrame(all_elements_expanded_gr, keep.extra.columns = T)

# original GRanges
all_elements_original_gr <- makeGRangesFromDataFrame(all_elements, keep.extra.columns = T)

################################################################################### filtering out elements which are in 150 kb of another element
# finding overlaps
all_elements_overlaps <- findOverlaps(all_elements_expanded_gr, all_elements_original_gr)

# combining both in data.frame, removing self-overlaps
all_elements_overlaps <- data.frame(fullName_original = all_elements_original_gr[subjectHits(all_elements_overlaps)]$fullName,
                                    fullName_expanded = all_elements_expanded_gr[queryHits(all_elements_overlaps)]$fullName, 
                                    stringsAsFactors = F)
all_elements_overlaps <- all_elements_overlaps[all_elements_overlaps$"fullName_original" != all_elements_overlaps$"fullName_expanded", ]
all_elements_overlaps <- unique(c(all_elements_overlaps$"fullName_expanded", all_elements_overlaps$"fullName_original"))

# filtering from original table 
all_elements_filtered <- all_elements[!(all_elements$"fullName" %in% all_elements_overlaps), ]

################################################################################### distance between elements and genes 
# finding distance between LTR elements and genes
all_elements_original_gr <- makeGRangesFromDataFrame(all_elements_filtered, keep.extra.columns = T)
all_elements_original_gr$distance_to_genes <- mcols(distanceToNearest(all_elements_original_gr, knownGenes, ignore.strand = T))$distance

##################################################################################### run from here
##################################################################################### 
# filtering based on: 
#   - element name (MT2_Mm, ORR1A, MT2_full)
#   - distance to nearest gene  
#   - hand picked filter 
#   - 100 elements ordered by FPKM in 2-cell stage 

# names picked by hand for filtering
hand_picked_filter_elements <- "chr5:151679590-151685986:+|MT2_full"

# hand_picked_filter_elements <- c("chr1:93965701-93966191:+|MT2_Mm", "chr2:157528079-157528349:+|MT2_Mm", "chr5:137721893-137722008:-|MT2_Mm",
#                                  "chr6:89202695-89202891:-|MT2_Mm", "chr12:19108326-19108819:-|MT2_Mm", "chr12:19108326-19108819:-|MT2_Mm",
#                                  "chr1:83031132-83031678:+|MT2_Mm", "chr17:6420849-6421203:+|ORR1A0", "chr17:6664850-6665198:-|ORR1A0",
#                                  "chr3:88577651-88577998:-|ORR1A0", "chr13:10040804-10041141:+|ORR1A0", "chr1:85172507-85178976:+|MT2_full",
#                                  "chr12:19894744-19901175:+|MT2_full", "chr13:76077210-76083662:-|MT2_full", "chr11:60651277-60657678:-|MT2_full",
#                                  "chr3:79046228-79052711:-|MT2_full", "chr12:19931732-19938183:-|MT2_full", "chr5:151679590-151685986:+|MT2_full")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/coverage_scanning/hand_filter")
all_element_names <- c("MT2_full", "MT2_Mm", "ORR1A0")

for(element_name in all_element_names){
  
  for(min_distance in c(0, 50000, 100000, 150000)){
    
    element_name <- "MT2_full"
    min_distance <- 0
    
    # original ranges
    element_original_ranges <- all_elements_original_gr[all_elements_original_gr$distance_to_genes >= min_distance & all_elements_original_gr$repName == element_name]
    element_original_ranges <- element_original_ranges[!(element_original_ranges$fullName %in% hand_picked_filter_elements)]
    element_original_ranges <- element_original_ranges[order(element_original_ranges$FPKM, decreasing = T)]
    
    # take top 100 based on FPKM in 2-cell
    if(length(element_original_ranges) >= 100){
      element_original_ranges <- element_original_ranges[1:100]
    }
    
    mcols(element_original_ranges)$FPKM <- NULL
    element_width <- max(width(element_original_ranges))
    
    # expand ranges 200kb up and downstream
    element_expanded_ranges <- as.data.frame(element_original_ranges)
    element_expanded_ranges$end <- element_expanded_ranges$start + 2e5 
    element_expanded_ranges$start <- element_expanded_ranges$start - 2e5 
    element_expanded_ranges <- makeGRangesFromDataFrame(element_expanded_ranges, keep.extra.columns = T)
    
    ##################################################################################### filtering with repeatMasker and/or knownGene
    # removing original element ranges from filter
    genomic_filter <- c(repeatMasker)
    genomic_filter <- GenomicRanges::setdiff(genomic_filter, element_original_ranges, ignore.strand = T)
    seqlevels(genomic_filter, force = T) <- seqlevels(bam_2C)
    
    # # setting file name
    # file_name <- paste0(element_name, "_repeatMaskerKnownGeneExonsFiltered")
    # file_name <- paste0(element_name, "_repeatMaskerKnownGenesFiltered")
    file_name <- paste0(element_name, "_repeatMaskerFiltered_", "geneDist", (min_distance / 1000), "kb_all")
    
    # filtering with genomic filter, getting data.frame and GRanges
    element_genomic_filtered <- lapply(element_expanded_ranges, function(x) GenomicRanges::setdiff(x, genomic_filter, ignore.strand = T))
    names(element_genomic_filtered) <- element_expanded_ranges$fullName
    element_genomic_filtered <- lapply(element_genomic_filtered, as.data.frame)
    element_genomic_filtered_df <- do.call(rbind, element_genomic_filtered)
    element_genomic_filtered_gr <- makeGRangesFromDataFrame(element_genomic_filtered_df)
    
    ##################################################################################### filtering with FPKM threshold 
    # getting counts of all samples with summarizeOverlaps function
    count_list <- lapply(list(bam_GV, bam_1C, bam_2C, bam_2C_aphi, bam_4C), countFunction, feature = element_genomic_filtered_gr)
    count_df <- as.data.frame(do.call(cbind, count_list))
    colnames(count_df) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
    count_df$feature_names <- rownames(count_df)
    rownames(count_df) <- NULL
    
    # calculating FPKM from counts
    element_genomic_filtered_fpkm <- cbind(element_genomic_filtered_df, count_df)
    element_genomic_filtered_fpkm <- element_genomic_filtered_fpkm[, -grep("strand", colnames(element_genomic_filtered_fpkm))]
    rownames(element_genomic_filtered_fpkm) <- NULL
    invisible(lapply(X = names(number_of_reads),
                     FUN = function(X) element_genomic_filtered_fpkm[, X] <<- element_genomic_filtered_fpkm[, X] / (number_of_reads[X] * (element_genomic_filtered_fpkm$width / 1000))))
    
    # removing MT2 positions from FPKM data.frame 
    element_genomic_filtered_fpkm_gr <- makeGRangesFromDataFrame(element_genomic_filtered_fpkm, keep.extra.columns = T)
    element_genomic_filtered_fpkm_gr <- element_genomic_filtered_fpkm_gr[-queryHits(findOverlaps(element_genomic_filtered_fpkm_gr, element_original_ranges))]
    element_genomic_filtered_fpkm <- element_genomic_filtered_fpkm[element_genomic_filtered_fpkm$feature_names %in% element_genomic_filtered_fpkm_gr$feature_names, ]
    
    # filtering ranges by FPKM values
    bam_names <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
    element_fpkm_filtered <- lapply(bam_names, 
                                    FUN = filterRangesByFPKM, 
                                    element_fpkm = element_genomic_filtered_fpkm, 
                                    element_ranges = element_genomic_filtered_df,
                                    fpkm_limit = 1)
    names(element_fpkm_filtered) <- c("bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
    
    #########################################################################################################
    # getting coverage as list of coverage vectors (with names = relative position) for each of 100 elements (5x, for each stage/.bam file)
    # GV
    element_coverage_GV <- lapply(X = 1:length(element_fpkm_filtered$bam_GV), 
                                  function(X) coverageForSummedPlot(X, coverage_stage = coverage_GV, 
                                                                    bam_stage_name = "bam_GV", 
                                                                    element_grList = element_fpkm_filtered)) 
    # 1-cell
    element_coverage_1C <- lapply(X = 1:length(element_fpkm_filtered$bam_1C), 
                                  function(X) coverageForSummedPlot(X, coverage_stage = coverage_1C, 
                                                                    bam_stage_name = "bam_1C", 
                                                                    element_grList = element_fpkm_filtered)) 
    # 2-cell
    element_coverage_2C <- lapply(X = 1:length(element_fpkm_filtered$bam_2C), 
                                  function(X) coverageForSummedPlot(X, coverage_stage = coverage_2C, 
                                                                    bam_stage_name = "bam_2C", 
                                                                    element_grList = element_fpkm_filtered)) 
    # 2-cell + aphidicolin
    element_coverage_2C_aphi <- lapply(X = 1:length(element_fpkm_filtered$bam_2C_aphi), 
                                       function(X) coverageForSummedPlot(X, coverage_stage = coverage_2C_aphi, 
                                                                         bam_stage_name = "bam_2C_aphi", 
                                                                         element_grList = element_fpkm_filtered)) 
    # 4-cell + aphidicolin
    element_coverage_4C <- lapply(X = 1:length(element_fpkm_filtered$bam_4C), 
                                  function(X) coverageForSummedPlot(X, coverage_stage = coverage_4C, 
                                                                    bam_stage_name = "bam_4C", 
                                                                    element_grList = element_fpkm_filtered)) 
    
    ######################################################################################################### summing coverage of all elements of one stage
    element_coverage_all <- lapply(list(element_coverage_GV, element_coverage_1C, element_coverage_2C, element_coverage_2C_aphi, element_coverage_4C), positionCoverageDF)
    element_coverage_all_df <- do.call(rbind, element_coverage_all)
    element_coverage_all_df$stage <- rep(c("GV", "1C", "2C", "2Ca", "4C"), each = (2 * 1.5e5 + element_width + 1))
    element_coverage_all_df$stage <- factor(element_coverage_all_df$stage, levels = c("GV", "1C", "2C", "2Ca", "4C"))
    element_coverage_all_df$count[element_coverage_all_df$count == 0] <- NA
    element_coverage_all_df$pos[is.na(element_coverage_all_df$count)] <- NA
    
    ######################################################################################################### plotting
    # coverage plot
    ggplot(element_coverage_all_df, aes(x = pos, y = count)) +
      geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
      scale_fill_manual(values = c("grey", "black"), guide = FALSE) +
      scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) + 
      coord_cartesian(ylim = c(0, 100)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      facet_grid(stage ~ .) +
      ggsave(filename = paste0(file_name, "_summedCoverage.pdf"), width = 30, height = 10)
    
    # bin plot
    element_coverage_merged <- lapply(X = 1:5, FUN = function(X) element_coverage_all[[X]][, 1:2])
    element_coverage_merged <- Reduce(function(...) merge(..., by = 'pos', all.x = TRUE), element_coverage_merged)
    colnames(element_coverage_merged) <- c("pos", "bam_GV", "bam_1C", "bam_2C", "bam_2C_aphi", "bam_4C")
    
    plotCummulativeFPKMBin(10000, merged_coverage_df = element_coverage_merged)
    plotCummulativeFPKMBin(2000, merged_coverage_df = element_coverage_merged)
    
  }
}
