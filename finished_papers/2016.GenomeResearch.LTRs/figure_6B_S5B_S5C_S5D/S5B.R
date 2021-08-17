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
library(tibble)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(BiocParallel)

######################################################## FUNCTIONS
countFunction <- function(bam_reads, input_features){
  
  register(MulticoreParam(workers = 8))
  se_filtered <- summarizeOverlaps(features = input_features,
                                   reads = bam_reads, 
                                   mode = "Union", 
                                   singleEnd = TRUE, 
                                   ignore.strand = TRUE)
  
  return(as.data.frame(assay(se_filtered)))
  
}

filterRangesByFPKM <- function(sample_name, fpkm_limit, fpkm_filtered_df, genomic_filtered_df){
  
  # get names of features above FPKM limit
  fpkm_above_limit_names <- fpkm_filtered_df[fpkm_filtered_df[, sample_name] > fpkm_limit, "feature_names"] 
  
  # set rownames to separate column 
  genomic_filtered_df %<>%
    tibble::rownames_to_column(var = "feature_name")
  
  # filter features above FPKM limit
  fpkm_filtered_grlist <- 
    genomic_filtered_df %>% 
    filter(!(feature_name %in% fpkm_above_limit_names)) %>% 
    mutate(feature_name = str_replace_all(feature_name, "\\.[0-9]*", "")) 
  
  # split to GRangesList
  fpkm_filtered_grlist %<>% 
    makeGRangesFromDataFrame(.) %>% 
    set_names(fpkm_filtered_grlist$feature_name) %>% 
    split(., names(.))
  
  return(fpkm_filtered_grlist)
  
}

coverageForSummedPlot <- function(x, sample_coverage, sample_name, element_grList){
  
  # function takes (filtered) ranges of one element and returns relative 
  # coverage on each position of expanded coordinates (+-2e5 kb) in .bam file
  
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
  names(counts_zero) <- start(element_expanded) : end(element_expanded)
  
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
  
  # wide format
  all_coverage_counts_df_sum_wide <- 
    merged_coverage_df %>% 
    dplyr::select(stage, pos, count, id) %>% 
    mutate(count = replace(count, is.na(count), 0)) %>% 
    tidyr::spread(key = stage, value = count) %>% 
    dplyr::select(-id) 
  
  # convert to binned data
  all_coverage_counts_df_sum_wide <- 
    rbind(all_coverage_counts_df_sum_wide %>% 
            dplyr::filter(pos >= -150000 & pos < 0) %>% 
            group_by(indx = gl(ceiling(nrow(.)/tile_width), tile_width, nrow(.))) %>%
            dplyr::select(-pos) %>%
            summarise_each(funs(sum)) %>% 
            mutate(position = paste(seq(150000, tile_width, -tile_width), "upstream", sep = "_")),
          all_coverage_counts_df_sum_wide %>% 
            dplyr::filter(pos >= 0 & pos <= element_width) %>% 
            dplyr::select(-pos) %>%
            summarise_each(funs(sum)) %>% 
            mutate(position = "element", 
                   indx = 1) %>% 
            dplyr::select(indx, s_Oo, s_2C, position),
          all_coverage_counts_df_sum_wide %>% 
            dplyr::filter(pos > element_width & pos <= (150000 + element_width)) %>% 
            group_by(indx = gl(ceiling(nrow(.)/tile_width), tile_width, nrow(.))) %>%
            dplyr::select(-pos) %>%
            summarise_each(funs(sum)) %>% 
            mutate(position = paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_"))) %>%
    dplyr::select(-1) %>% 
    mutate(s_Oo = s_Oo / (number_of_reads["s_Oo"] * (tile_width / 1000)), 
           s_2C = s_2C / (number_of_reads["s_2C"] * (tile_width / 1000)))
  
  # output table
  all_coverage_counts_df_sum_wide %>% 
    mutate(pos = c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width))) %>% 
    write_csv(path = paste0("Park2013_", choosen_element_name, "_cummulativeFPKM_bins", tile_width, ".csv"), col_names = T)
  
  # long format
  all_coverage_counts_df_sum_long <- 
    all_coverage_counts_df_sum_wide %>% 
    tidyr::gather(key = stage, value = fpkm, -position) %>% 
    mutate(element = ifelse(position == "element", "in_element",  "out_element"), 
           width = ifelse(position == "element", element_width, tile_width), 
           pos = rep(c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width)), length(number_of_reads)), 
           stage = factor(stage, levels = names(number_of_reads))) %>% 
    dplyr::select(pos, fpkm, element, stage, position, width) %>% 
    arrange(desc(stage))
  
  # bin plot
  ggplot(all_coverage_counts_df_sum_long, aes(x = pos, y = fpkm)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + width, ymax = fpkm, fill = stage)) +
    scale_fill_manual(values = c("black", "grey30")) +
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
    ggsave(filename = paste0("Park2013_", choosen_element_name, "_cummulativeFPKM_bins", tile_width, "_binPlot.pdf"), width = 30, height = 10)
  
}

######################################################## READ DATA
# .bam files
bam_Oo <- readGAlignments("s_Oo.SE.bam")
bam_2C <- readGAlignments("s_1C.SE.bam")

# coverage
coverage_Oo <- coverage(bam_Oo)
coverage_2C <- coverage(bam_2C)

# repeatMasker (UCSC)
repeatMasker <- 
  read_delim("UCSC_repeatMasker_mm9_20170205_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, gene_id = repName) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# knownGene (UCSC)
knownGenes <- 
  makeTxDbFromGFF("UCSC_knownGene_mm9_20170205.gtf.gz") %>% 
  genes(.)

# combined repeatMasker and knownGene
genomic_filter <- c(knownGenes, repeatMasker)

# library size
number_of_reads <- 
  read_lines("Park_2013_library_size.txt") %>% 
  as.integer() %>% 
  divide_by(10^6) %>% 
  set_names(c("s_Oo", "s_2C"))

### coordinates and FPKM expression in 2-cell stage 
# MT2/ORR1A0 solo LTRs 
MT2_ORR1A0_all <- 
  read_delim("MT2_ORR1A0_solo_LTRs_orderedByExpressionIn2cell_Park_mm9.txt", delim = "\t") %>% 
  dplyr::rename(element_name = repName)

# MuERV-L elements
MT2_full <- 
  read_delim("MT2_full_orderedByExpressionIn2cell_Park_mm9.txt", delim = "\t") %>% 
  dplyr::select(-id) %>% 
  mutate(element_name = "MT2_full")

######################################################## MAIN CODE
### all elements ranges and expanded ranges
# combined elements
all_elements <- 
  rbind(MT2_ORR1A0_all, MT2_full) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name)) %>% 
  dplyr::filter(!(fullName %in% c("chr18:90631939-90632457:+|MT2_Mm"))) # filters out elements which go out of range of chromosome when expanded

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
all_elements_overlaping_names <- 
  cbind(all_elements_original_gr[subjectHits(all_elements_overlaps)]$fullName, 
        all_elements_expanded_gr[queryHits(all_elements_overlaps)]$fullName) %>% 
  as.data.frame(.) %>% 
  set_colnames(c("fullName_original", "fullName_expanded")) %>% 
  dplyr::filter(fullName_original != fullName_expanded) %$% 
  fullName_expanded %>% 
  unique() %>% 
  as.character()

# filter from original table  
all_elements %<>% 
  dplyr::filter(!(fullName %in% all_elements_overlaping_names))


### summed coverage of top 100 elements, plot
# top 100 by expression for each element class
all_elements_top100 <- 
  all_elements %>%
  group_by(element_name) %>%
  arrange(desc(s_2C_fpkm)) %>%
  dplyr::slice(1:100) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-s_2C_fpkm) 

# loop through 3 element types
for(choosen_element_name in unique(all_elements_top100$element_name)){
  
  ### getting elements of one type
  # original element coordinates
  element_original <- 
    all_elements_top100 %>% 
    dplyr::filter(element_name == choosen_element_name) %>% 
    dplyr::select(-width)
  element_original_ranges <- makeGRangesFromDataFrame(element_original, keep.extra.columns = T)
  element_width <- max(width(element_original_ranges))
  
  # expand coordinates 200kb up and downstream
  element_expanded <- 
    element_original %>% 
    mutate(end = start + 2e5, 
           start = start - 2e5) 
  element_expanded_ranges <- makeGRangesFromDataFrame(element_expanded, keep.extra.columns = T)
  
  
  ### filter with repeatMasker and knownGene UCSC tables
  # remove coordinates of elements from filter table
  genomic_filter_without_elements <- GenomicRanges::setdiff(genomic_filter, element_original_ranges, ignore.strand = T)
  seqlevels(genomic_filter_without_elements, force = T) <- seqlevels(bam_2C)
  
  # filter
  element_genomic_filtered_df <- 
    lapply(element_expanded_ranges, function(x) GenomicRanges::setdiff(x, genomic_filter_without_elements, ignore.strand = T)) %>% 
    set_names(element_original$fullName) %>% 
    lapply(., as.data.frame) %>% 
    do.call(rbind, .) 
  element_genomic_filtered_ranges <- makeGRangesFromDataFrame(element_genomic_filtered_df)
  
  
  ### filter with FPKM threshold 
  # counts all samples
  count_list <- lapply(list(bam_Oo, bam_2C), countFunction, input_features = element_genomic_filtered_ranges)
  
  # FPKM
  element_fpkm_df <- 
    do.call(cbind, count_list) %>% 
    set_colnames(c("s_Oo", "s_2C")) %>% 
    mutate(feature_names = rownames(.)) %>% 
    set_rownames(NULL) %>% 
    cbind(element_genomic_filtered_df, .) %>% 
    set_rownames(NULL) %>% 
    mutate(s_Oo = (s_Oo / (number_of_reads["s_Oo"] * (width / 1000))), 
           s_2C = (s_2C / (number_of_reads["s_2C"] * (width / 1000))))
  
  # remove coordinates of elements from FPKM table 
  element_fpkm_ranges <- makeGRangesFromDataFrame(element_fpkm_df, keep.extra.columns = T)
  element_fpkm_ranges <- element_fpkm_ranges[-queryHits(findOverlaps(element_fpkm_ranges, element_original_ranges))]
  element_fpkm_filtered_df <- element_fpkm_df[element_fpkm_df$feature_names %in% element_fpkm_ranges$feature_names, ]
  
  # filter
  element_fpkm_filtered_grList <- 
    lapply(names(number_of_reads), 
           FUN = filterRangesByFPKM, 
           fpkm_limit = 1, 
           fpkm_filtered_df = element_fpkm_filtered_df, 
           genomic_filtered_df = element_genomic_filtered_df) %>% 
    set_names(names(number_of_reads))
  
  
  ### coverage of expanded top 100 elements in each .bam file
  # oocyte
  element_coverage_Oo <- lapply(X = 1:length(element_fpkm_filtered_grList$s_Oo), 
                                function(X) coverageForSummedPlot(X, 
                                                                  sample_coverage = coverage_Oo, 
                                                                  sample_name = "s_Oo", 
                                                                  element_grList = element_fpkm_filtered_grList)) 
  # 2-cell
  element_coverage_2C <- lapply(X = 1:length(element_fpkm_filtered_grList$s_2C), 
                                function(X) coverageForSummedPlot(X, 
                                                                  sample_coverage = coverage_2C, 
                                                                  sample_name = "s_2C", 
                                                                  element_grList = element_fpkm_filtered_grList)) 
  
  
  ### sum coverage, output table, plot coverage and binned coverage
  # sum coverage of all elements in one .bam
  element_coverage_all_df <- 
    lapply(list(element_coverage_Oo, element_coverage_2C), positionCoverageDF) %>% 
    bind_rows(.) %>% 
    mutate(stage = rep(c("s_Oo", "s_2C"), each = (2 * 1.5e5 + element_width + 1)), 
           stage = factor(stage, levels = c("s_Oo", "s_2C")), 
           count = replace(count, count == 0, NA), 
           id = rep(1:(nrow(.)/2), length(number_of_reads)))
  
  # plot coverage
  ggplot(element_coverage_all_df, aes(x = pos, y = count)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    scale_fill_manual(values = c("grey", "black")) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_grid(stage ~ .) +
    ggsave(filename = paste0("Park2013_", choosen_element_name, "_summedCoverage.pdf"), width = 30, height = 10)
  
  # bin plot
  plotCummulativeFPKMBin(merged_coverage_df = element_coverage_all_df, tile_width = 10000)
  
}