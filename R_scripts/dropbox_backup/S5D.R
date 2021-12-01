rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
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

######################################################## FUNCTIONS
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

coverageForSummedPlot <- function(x, sample_coverage, sample_name, element_grList, element_original_ranges_fun, element_expanded_ranges_fun){
  
  # function takes (filtered) ranges of one element and returns relative 
  # coverage on each position of expanded range (+-2e5 kb) in .bam file
  
  # take element ranges, get original element data
  element_gr <- element_grList[[sample_name]][[x]]
  element_original <- element_original_ranges_fun[element_original_ranges_fun$fullName == names(element_gr)[1]]
  element_expanded <- element_expanded_ranges_fun[element_original_ranges_fun$fullName == names(element_gr)[1]]
  
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

positionCoverageDF <- function(element_coverage, element_width_fun){
  
  # get all unique position across all elements in one stage
  all_positions <- unique(unlist(lapply(element_coverage, names)))
  
  # position matrix for all counts (rows = position, columns = elements)
  position_matrix <- matrix(NA, 
                            nrow = length(all_positions), 
                            ncol = length(element_coverage), 
                            dimnames = list(all_positions, c(1:length(element_coverage))))
  
  # fill position matrix with coverage based on position names
  for (i in seq_along(element_coverage)) {
    position_matrix[names(element_coverage[[i]]), i] <- element_coverage[[i]]
  }
  
  # keep only position for which all elements have position and then keep only +-1.5e5 kb
  position_matrix <- position_matrix[complete.cases(position_matrix), ]
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000") : which(rownames(position_matrix) == as.character(150000 + element_width_fun)), ]
  
  # sum rows = summed coverage of all elements per position
  element_coverage_summed <- 
    rowSums(position_matrix) %>% 
    data.frame(pos = as.numeric(names(.)), 
               count = .) %>% 
    mutate(element = ifelse(pos >= 0 & pos <= element_width_fun, "in_element", "out_element"))
  
  return(element_coverage_summed)
  
}

######################################################## READ DATA
# .bam files
bam_GV <- readGAlignmentPairs("s_GV.WE.bam")
bam_1C <- readGAlignmentPairs("s_1cell.WE.bam")
bam_2C <- readGAlignmentPairs("s_2cell.WE.bam")
bam_2Ca <- readGAlignmentPairs("s_2cell.WE_DNAm.bam")

# coverage
coverage_GV <- coverage(bam_GV)
coverage_1C <- coverage(bam_1C)
coverage_2C <- coverage(bam_2C)
coverage_2Ca <- coverage(bam_2Ca)

# repeatMasker (UCSC)
repeatMasker <- 
  read_delim("UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, element_name = repName) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# knownGene (UCSC)
knownGenesTxDb <- makeTxDbFromGFF("UCSC_knownGene_mm10_20161126.gtf.gz")
knownGenes <- genes(knownGenesTxDb)
names(mcols(knownGenes)) <- "element_name"

# library size (from STAR aligner Log.final.out)
logs_filenames <- file.path(c("s_GV.WELog.final.out", "s_1cell.WELog.final.out", "s_2cell.WELog.final.out", "s_2cell.WE_DNAmLog.final.out"))
library_size_df <- 
  sapply(X = 1:length(logs_filenames), 
         FUN = function(X){
           read_delim(logs_filenames[X], delim = "\t", col_names = F)[8, 2] %>%
             as.integer(.)
         }) %>% 
  data.frame(library_size = . / 10^6) %>% 
  mutate(stage = factor(c("GV", "1C", "2C", "2Ca"), levels = c("GV", "1C", "2C", "2Ca"))) %>% 
  dplyr::select(stage, library_size) 

### coordinates and FPKM expression in 2-cell stage 
# MT2/ORR1A0 solo LTRs 
MT2_ORR1A0_all <- 
  read_delim("MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", delim = "\t") %>% 
  dplyr::select(-4)

# MuERV-L elements 
MT2_full <- 
  read_delim("MT2_full_orderedByFPKMIn2cell.txt", delim = "\t") %>% 
  dplyr::select(-c(4, 8)) %>% 
  mutate(repName = "MT2_full")

######################################################## MAIN CODE
### all elements ranges and expanded ranges
# combined elements
all_elements <- 
  rbind(MT2_ORR1A0_all, MT2_full) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           repName)) %>% 
  dplyr::filter(!grepl("random|Un", seqnames), 
                !(fullName %in% c("chr5:151679590-151685986:+|MT2_full"))) # filter out elements which go out of range of chromosome when expanded

# original coordinates
all_elements_original_gr <- makeGRangesFromDataFrame(all_elements, keep.extra.columns = T)

# expand ranges by 150kb in both directions
all_elements_expanded_gr <- 
  all_elements %>% 
  mutate(end = start + 1.5e5, 
         start = start - 1.5e5) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)


### filter out elements which are closer than 150 kb to another element
# find overlaps
all_elements_overlaps <- findOverlaps(all_elements_expanded_gr, all_elements_original_gr)

# remove self-overlaps
all_elements_overlaps <- data.frame(fullName_original = all_elements_original_gr[subjectHits(all_elements_overlaps)]$fullName,
                                    fullName_expanded = all_elements_expanded_gr[queryHits(all_elements_overlaps)]$fullName, 
                                    stringsAsFactors = F)
all_elements_overlaps <- all_elements_overlaps[all_elements_overlaps$"fullName_original" != all_elements_overlaps$"fullName_expanded", ]
all_elements_overlaps <- unique(c(all_elements_overlaps$"fullName_expanded", all_elements_overlaps$"fullName_original"))

# filter from original table 
all_elements_original_gr <- 
  all_elements[!(all_elements$"fullName" %in% all_elements_overlaps), ] %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# distance between elements and genes
all_elements_original_gr$distance_to_genes <- mcols(distanceToNearest(all_elements_original_gr, knownGenes, ignore.strand = T))$distance


### summed coverage of top 100 elements, plot
# loop through 3 element types
element_coverage_raw_list <- lapply(c("MT2_full", "MT2_Mm", "ORR1A0"), function(element_name){
  
  # loop through minimal distance to closest gene
  lapply(c(0, 50000), function(min_distance){
    
    # set file name
    file_name <- paste0(element_name, "_repeatMaskerFiltered_", min_distance / 1000, "kb")
    
    ### top 100 elements of one type, original and expanded coordinates
    # original element coordinates
    element_original_ranges <- all_elements_original_gr[all_elements_original_gr$distance_to_genes >= min_distance & all_elements_original_gr$repName == element_name]
    
    # top 100 by FPKM in 2-cell
    element_original_ranges <- element_original_ranges[order(element_original_ranges$FPKM, decreasing = T)]
    element_original_ranges <- element_original_ranges[1:100]
    mcols(element_original_ranges)$FPKM <- NULL
    element_width <- max(width(element_original_ranges))
    
    # expand ranges 200kb up and downstream
    element_expanded_ranges <- 
      as.data.frame(element_original_ranges) %>% 
      dplyr::mutate(end = start + 2e5, 
                    start = start - 2e5) %>% 
      makeGRangesFromDataFrame(., keep.extra.columns = T)
    
    
    ### filter with repeatMasker and knownGene UCSC tables
    # remove original element ranges from filter table
    genomic_filter <- c(repeatMasker)
    genomic_filter <- GenomicRanges::setdiff(genomic_filter, element_original_ranges, ignore.strand = T)
    seqlevels(genomic_filter, force = T) <- seqlevels(bam_2C)
    
    # filter
    element_genomic_filtered_df <- 
      lapply(element_expanded_ranges, function(x) GenomicRanges::setdiff(x, genomic_filter, ignore.strand = T)) %>% 
      set_names(element_expanded_ranges$fullName) %>% 
      lapply(., as.data.frame) %>% 
      do.call(rbind, .) %>% 
      mutate(feature_name = rownames(.)) %>% 
      set_rownames(NULL)
    
    
    ### filter with FPKM threshold 
    # counts in all samples
    element_genomic_filtered_fpkm <- 
      lapply(list(bam_GV, bam_1C, bam_2C, bam_2Ca), countFunction, feature = makeGRangesFromDataFrame(element_genomic_filtered_df)) %>% 
      bind_cols(.) %>% 
      set_colnames(c("GV", "1C", "2C", "2Ca")) %>% 
      cbind(element_genomic_filtered_df, .) %>% 
      dplyr::select(-strand) 
    
    # FPKM
    invisible(lapply(X = as.character(library_size_df$stage), 
                     FUN = function(X){
                       element_genomic_filtered_fpkm[, X] <<- 
                         element_genomic_filtered_fpkm[, X] / (library_size_df[library_size_df$stage == X, "library_size"] * (element_genomic_filtered_fpkm$width / 1000))
                     }))
    
    # remove element positions from FPKM table 
    element_genomic_filtered_fpkm_gr <- makeGRangesFromDataFrame(element_genomic_filtered_fpkm, keep.extra.columns = T)
    element_genomic_filtered_fpkm_gr <- element_genomic_filtered_fpkm_gr[-queryHits(findOverlaps(element_genomic_filtered_fpkm_gr, element_original_ranges))]
    element_genomic_filtered_fpkm <- element_genomic_filtered_fpkm[element_genomic_filtered_fpkm$feature_name %in% element_genomic_filtered_fpkm_gr$feature_name, ]
    colnames(element_genomic_filtered_fpkm)[6:9] <- paste0("bam_", colnames(element_genomic_filtered_fpkm)[6:9])
    
    # filter
    element_fpkm_filtered <-
      lapply(c("bam_GV", "bam_1C", "bam_2C", "bam_2Ca"), 
             FUN = filterRangesByFPKM, 
             element_fpkm = element_genomic_filtered_fpkm, 
             element_ranges = element_genomic_filtered_df,
             fpkm_limit = 1) %>% 
      set_names(c("bam_GV", "bam_1C", "bam_2C", "bam_2Ca"))
    
    
    ### coverage of expanded top 100 elements in each .bam file
    # GV
    element_coverage_GV <- lapply(X = 1:length(element_fpkm_filtered$bam_GV), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_GV, 
                                                                    sample_name = "bam_GV", 
                                                                    element_grList = element_fpkm_filtered, 
                                                                    element_original_ranges_fun = element_original_ranges, 
                                                                    element_expanded_ranges_fun = element_expanded_ranges)) 
    # 1-cell
    element_coverage_1C <- lapply(X = 1:length(element_fpkm_filtered$bam_1C), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_1C, 
                                                                    sample_name = "bam_1C", 
                                                                    element_grList = element_fpkm_filtered, 
                                                                    element_original_ranges_fun = element_original_ranges, 
                                                                    element_expanded_ranges_fun = element_expanded_ranges)) 
    # 2-cell
    element_coverage_2C <- lapply(X = 1:length(element_fpkm_filtered$bam_2C), 
                                  function(X) coverageForSummedPlot(X, sample_coverage = coverage_2C, 
                                                                    sample_name = "bam_2C", 
                                                                    element_grList = element_fpkm_filtered, 
                                                                    element_original_ranges_fun = element_original_ranges, 
                                                                    element_expanded_ranges_fun = element_expanded_ranges)) 
    # 2-cell + aphidicolin
    element_coverage_2Ca <- lapply(X = 1:length(element_fpkm_filtered$bam_2Ca), 
                                   function(X) coverageForSummedPlot(X, sample_coverage = coverage_2Ca, 
                                                                     sample_name = "bam_2Ca", 
                                                                     element_grList = element_fpkm_filtered, 
                                                                     element_original_ranges_fun = element_original_ranges, 
                                                                     element_expanded_ranges_fun = element_expanded_ranges)) 
    
    
    ### sum coverage, plot
    # sum coverage of all elements in one .bam
    element_coverage_all_df <- 
      lapply(list(element_coverage_GV, element_coverage_1C, element_coverage_2C, element_coverage_2Ca), positionCoverageDF, element_width_fun = element_width) %>% 
      bind_rows(.) %>% 
      mutate(stage = rep(c("GV", "1C", "2C", "2Ca"), each = (2 * 1.5e5 + element_width + 1)), 
             stage = factor(stage, levels = c("GV", "1C", "2C", "2Ca")), 
             count = replace(count, count == 0, NA), 
             pos = replace(pos, is.na(count), NA), 
             element_name = element_name, 
             min_distace = min_distance, 
             element_width = element_width)
    
    # # coverage plot
    # ggplot(element_coverage_all_df, aes(x = pos, y = count)) +
    #   geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    #   scale_fill_manual(values = c("grey", "black"), guide = FALSE) +
    #   scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) + 
    #   coord_cartesian(ylim = c(0, 100)) +
    #   theme_bw() +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank()) +
    #   facet_grid(stage ~ .) +
    #   ggsave(filename = paste0(file_name, "_summedCoverage.pdf"), width = 30, height = 10)
    
    return(element_coverage_all_df)
    
  })
  
})

# unlist, save as .RDS
element_coverage_raw_list_unlisted <- unlist(element_coverage_raw_list, recursive = FALSE)
saveRDS(object = element_coverage_raw_list_unlisted, file = "figure_S5D_coverage_list.RDS")
