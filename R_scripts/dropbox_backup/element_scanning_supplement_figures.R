library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(BiocParallel)

# options(bitmapType = "cairo")
# rm(list = ls()[!grepl("^bam_|^coverage_|^repeatMasker$|^knownGenes|countFunction|filterRangesByFPKM|coverageForSummedPlot|positionCoverageDF|plotCummulativeFPKMBin|library_size_df", ls())])

################################################################################### functions
countFunction <- function(feature, bam_reads){
  
  register(MulticoreParam(workers = 8))
  se_filtered <- summarizeOverlaps(features = feature,
                                   reads = bam_reads, 
                                   mode = "Union", 
                                   singleEnd = FALSE, 
                                   ignore.strand = TRUE)
  
  return(as.data.frame(assay(se_filtered)))
  
}

filterRangesByFPKM <- function(bam_name, element_fpkm, element_ranges, fpkm_limit){
  
  # get names of features above FPKM limit
  fpkm_df_names_limit <- element_fpkm[element_fpkm[, bam_name] > fpkm_limit, "feature_name"]
  
  # filter features above FPKM limit, split to GRangesList
  element_ranges_fpkm_filtered_grList <- 
    element_ranges %>% 
    dplyr::filter(!(feature_name %in% fpkm_df_names_limit)) %>% 
    mutate(feature_name = str_replace_all(feature_name, "\\.[0-9]*", "")) 
  
  element_ranges_fpkm_filtered_grList %<>% 
    makeGRangesFromDataFrame(.) %>% 
    set_names(element_ranges_fpkm_filtered_grList$feature_name) %>% 
    split(., names(.))
  
  return(element_ranges_fpkm_filtered_grList)
  
}

coverageForSummedPlot <- function(x, sample_coverage, sample_name, element_grList){
  
  # function takes (filtered) ranges of one element and returns relative 
  # coverage on each position of expanded range (+-2e5 kb) in .bam file
  
  # take element ranges, get original element data
  element_gr <- element_grList[[sample_name]][[x]]
  element_original <- element_original_ranges[element_original_ranges$fullName == names(element_gr)[1]]
  element_expanded <- element_expanded_ranges[element_original_ranges$fullName == names(element_gr)[1]]
  
  # find coverage
  element_coverage <- Views(sample_coverage, as(element_gr, "RangesList"))
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
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000") : which(rownames(position_matrix) == as.character(150000 + element_width)), ]
  
  # sum rows = summed coverage of all elements per position
  element_coverage_summed <- 
    rowSums(position_matrix) %>% 
    data.frame(pos = as.numeric(names(.)), 
               count = .) %>% 
    mutate(element = ifelse(pos >= 0 & pos <= element_width, "in_element", "out_element"))

  return(element_coverage_summed)
  
}

plotCummulativeFPKMBin <- function(tile_width, merged_coverage_df){
  
  # wide format using tidyr::spread
  # syntax:
  # - key = values to put as new columns
  # - value = values to fill cells 
  all_coverage_counts_df_sum_wide <- 
    merged_coverage_df %>% 
    dplyr::select(stage, pos, count) %>% 
    mutate(count = replace(count, is.na(count), 0)) %>% 
    tidyr::spread(key = stage, value = count) %>% 
    set_colnames(c("pos", "bam_GV", "bam_1C", "bam_2C", "bam_2Ca", "bam_4C"))
  
  all_coverage_counts_df_sum_wide <- 
    rbind(all_coverage_counts_df_sum_wide %>% #upstream
            filter(pos >= -150000 & pos < 0) %>%
            group_by(indx = gl(ceiling(n()/tile_width), tile_width, n())) %>%
            dplyr::select(2:6) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = paste(seq(150000, tile_width, -tile_width), "upstream", sep = "_")),
          all_coverage_counts_df_sum_wide %>% #element
            filter(pos >= 0 & pos <= element_width) %>%
            select(2:6) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = "element", 
                   indx = 1) %>%
            dplyr::select(c(7, 1:6)), 
          all_coverage_counts_df_sum_wide %>% #downstream
            filter(pos > element_width & pos <= (150000 + element_width)) %>%
            group_by(indx = gl(ceiling(n()/tile_width), tile_width, n())) %>%
            dplyr::select(2:6) %>%
            summarise_each(funs(sum)) %>%
            mutate(position = paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_"))) %>%
    select(-1)
  
  # normalizing  FPM to FPKM (for bins)
  all_coverage_counts_df_sum_wide[, 1:5] <- all_coverage_counts_df_sum_wide[, 1:5] / (tile_width / 1000)
  
  # output .csv
  all_coverage_counts_df_sum_wide %>%
    mutate(pos = c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width))) %>%
    write.csv(file = paste0(file_name, "_cummulativeFPKM_bins", tile_width, ".csv"), row.names = F)
  
  # long format for bin plot using tidyr::gather
  # syntax:
  # - gather all data except position (-position)
  # - name new key column stage (key = stage)
  # - name new value column fpkm (value = fpkm)
  all_coverage_counts_df_sum_long <- 
    all_coverage_counts_df_sum %>% 
    tidyr::gather(key = stage, value = fpkm, -position) %>% 
    mutate(element = ifelse((position == "element"), "in_element",  "out_element"), 
           width = ifelse((position == "element"), element_width, tile_width), 
           pos = rep(c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width)), 5), 
           stage = str_replace_all(stage, "bam_", "")) %>%
    select(pos, fpkm, element, stage, position, width) %>%
    mutate(stage = factor(stage, levels = c("2Ca", "2C", "1C", "4C", "GV"))) %>%
    arrange(stage) %>%
    mutate(stage = factor(stage, levels = c("GV", "4C", "1C", "2C", "2Ca")))
  
  # bin plot
  ggplot(all_coverage_counts_df_sum_long, aes(x = pos, y = fpkm)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + width, ymax = fpkm, fill = stage)) +
    scale_fill_manual(values = c("black", "grey60", "grey30", "grey50", "orange"), 
                      breaks = c("GV", "1C", "2C", "2Ca", "4C")) +
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
# setting working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/new/")

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
bam_2Ca <- readGAlignmentPairs(filenames[4])
bam_4C <- readGAlignmentPairs(filenames[5])

# coverage of BAM files
coverage_GV <- coverage(bam_GV)
coverage_1C <- coverage(bam_1C)
coverage_2C <- coverage(bam_2C)
coverage_2Ca <- coverage(bam_2Ca)
coverage_4C <- coverage(bam_4C)

# repeatMasker table
repeatMasker <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, gene_id = repName) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

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

library_size_df <- 
  sapply(X = 1:length(logs_filenames), 
         FUN = function(X){
           read_delim(logs_filenames[X], delim = "\t", col_names = F)[8, 2] %>%
             as.integer(.)
         }) %>% 
  data.frame(library_size = . / 10^6) %>% 
  mutate(stage = factor(c("GV", "1C", "2C", "2Ca", "4C"), levels = c("GV", "1C", "2C", "2Ca", "4C"))) %>% 
  dplyr::select(stage, library_size) 
       
# reading original element lists (MT2 and ORR1A0 solo LTRs ordered by expression in 2Cell)
MT2_ORR1A0_all <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", delim = "\t") %>% 
  dplyr::select(-4)

MT2_full <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_full_orderedByFPKMIn2cell.txt", delim = "\t") %>% 
  dplyr::select(-c(4, 8)) %>% 
  mutate(repName = "MT2_full")

################################################################################### all elements ranges and expanded ranges
# table with all elements (MT2/ORR1A solo LTRs, full MT2 elements)
all_elements <- 
  rbind(MT2_ORR1A0_all, MT2_full) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           repName)) %>% 
  dplyr::filter(!grepl("random|Un", seqnames))

# expanding ranges by 150kb in both directions, making expanded GRanges
all_elements_expanded_gr <- 
  all_elements %>% 
  mutate(end = start + 1.5e5, 
         start = start - 1.5e5) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

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

for(element_name in c("MT2_full", "MT2_Mm", "ORR1A0")){
  
  for(min_distance in c(0, 50000, 100000, 150000)){
    
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
    
    # setting file name
    # file_name <- paste0(element_name, "_repeatMaskerKnownGeneExonsFiltered")
    # file_name <- paste0(element_name, "_repeatMaskerKnownGenesFiltered")
    file_name <- paste0(element_name, "_repeatMaskerFiltered_", "geneDist", (min_distance / 1000), "kb_all")
    
    # filtering with genomic filter, getting data.frame and GRanges
    element_genomic_filtered_df <- 
      lapply(element_expanded_ranges, function(x) GenomicRanges::setdiff(x, genomic_filter, ignore.strand = T)) %>% 
      set_names(element_expanded_ranges$fullName) %>% 
      lapply(., as.data.frame) %>% 
      do.call(rbind, .) %>% 
      mutate(feature_name = rownames(.)) %>% 
      set_rownames(NULL)
    
    ##################################################################################### filtering with FPKM threshold 
    # getting counts of all samples with summarizeOverlaps
    element_genomic_filtered_fpkm <- 
      lapply(list(bam_GV, bam_1C, bam_2C, bam_2Ca, bam_4C), countFunction, feature = makeGRangesFromDataFrame(element_genomic_filtered_df)) %>% 
      bind_cols(.) %>% 
      set_colnames(c("GV", "1C", "2C", "2Ca", "4C")) %>% 
      cbind(element_genomic_filtered_df, .) %>% 
      dplyr::select(-strand) 
      
    # calculating FPKM from counts
    invisible(lapply(X = as.character(library_size_df$stage), 
                     FUN = function(X){
                       element_genomic_filtered_fpkm[, X] <<- 
                         element_genomic_filtered_fpkm[, X] / (library_size_df[library_size_df$stage == X, "library_size"] * (element_genomic_filtered_fpkm$width / 1000))
                     }))
    
    # removing MT2 positions from FPKM data.frame 
    element_genomic_filtered_fpkm_gr <- makeGRangesFromDataFrame(element_genomic_filtered_fpkm, keep.extra.columns = T)
    element_genomic_filtered_fpkm_gr <- element_genomic_filtered_fpkm_gr[-queryHits(findOverlaps(element_genomic_filtered_fpkm_gr, element_original_ranges))]
    element_genomic_filtered_fpkm <- element_genomic_filtered_fpkm[element_genomic_filtered_fpkm$feature_name %in% element_genomic_filtered_fpkm_gr$feature_name, ]
    colnames(element_genomic_filtered_fpkm)[6:10] <- paste0("bam_", colnames(element_genomic_filtered_fpkm)[6:10])
    
    # filtering ranges by FPKM values
    element_fpkm_filtered <-
      lapply(c("bam_GV", "bam_1C", "bam_2C", "bam_2Ca", "bam_4C"), 
             FUN = filterRangesByFPKM, 
             element_fpkm = element_genomic_filtered_fpkm, 
             element_ranges = element_genomic_filtered_df,
             fpkm_limit = 1) %>% 
      set_names(c("bam_GV", "bam_1C", "bam_2C", "bam_2Ca", "bam_4C"))

    #########################################################################################################
    # getting coverage as list of coverage vectors (with names = relative position) for each of 100 elements (5x, for each stage/.bam file)
    # GV
    element_coverage_GV <- lapply(X = 1:length(element_fpkm_filtered$bam_GV), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_GV, 
                                                                    sample_name = "bam_GV", 
                                                                    element_grList = element_fpkm_filtered)) 
    # 1-cell
    element_coverage_1C <- lapply(X = 1:length(element_fpkm_filtered$bam_1C), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_1C, 
                                                                    sample_name = "bam_1C", 
                                                                    element_grList = element_fpkm_filtered)) 
    # 2-cell
    element_coverage_2C <- lapply(X = 1:length(element_fpkm_filtered$bam_2C), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_2C, 
                                                                    sample_name = "bam_2C", 
                                                                    element_grList = element_fpkm_filtered)) 
    # 2-cell + aphidicolin
    element_coverage_2Ca <- lapply(X = 1:length(element_fpkm_filtered$bam_2Ca), 
                                       function(X) coverageForSummedPlot(X, sample_coverage = coverage_2Ca, 
                                                                         sample_name = "bam_2Ca", 
                                                                         element_grList = element_fpkm_filtered)) 

    # 4-cell + aphidicolin
    element_coverage_4C <- lapply(X = 1:length(element_fpkm_filtered$bam_4C), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_4C, 
                                                                    sample_name = "bam_4C", 
                                                                    element_grList = element_fpkm_filtered)) 
    
    ######################################################################################################### 
    # summing coverage of all elements of one stage
    element_coverage_all_df <- 
      lapply(list(element_coverage_GV, element_coverage_1C, element_coverage_2C, element_coverage_2Ca, element_coverage_4C), positionCoverageDF) %>% 
      bind_rows(.) %>% 
      mutate(stage = rep(c("GV", "1C", "2C", "2Ca", "4C"), each = (2 * 1.5e5 + element_width + 1)), 
             stage = factor(stage, levels = c("GV", "1C", "2C", "2Ca", "4C")), 
             count = replace(count, count == 0, NA), 
             pos = replace(pos, is.na(count), NA))

    # normalize for library size
    element_coverage_all_df_fpm <- 
      element_coverage_all_df %>% 
      left_join(., library_size_df, by = "stage") %>% 
      mutate(count = count / library_size) %>% 
      dplyr::select(-library_size)

    ######################################################################################################### plotting
    # coverage plot for raw data
    ggplot(element_coverage_all_df, aes(x = pos, y = count)) +
      geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
      scale_fill_manual(values = c("grey", "black"), guide = FALSE) +
      scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) + 
      coord_cartesian(ylim = c(0, 25)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      facet_grid(stage ~ .) +
      ggsave(filename = paste0(file_name, "_summedCoverage_raw_lim25.pdf"), width = 30, height = 10)
    
    # coverage plot for normalized data
    ggplot(element_coverage_all_df_fpm, aes(x = pos, y = count)) +
      geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
      scale_fill_manual(values = c("grey", "black"), guide = FALSE) +
      scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) + 
      coord_cartesian(ylim = c(0, 5)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      facet_grid(stage ~ .) +
      ggsave(filename = paste0(file_name, "_summedCoverage_fpm.pdf"), width = 30, height = 10)
    
    # bin plot and .csv output
    plotCummulativeFPKMBin(10000, merged_coverage_df = element_coverage_all_df_fpm)

    cat(file_name, "done", "\n")
    
  }
}
