library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)
library(broom)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(BiocParallel)

# options(bitmapType = 'cairo')
# rm(list = ls()[!grepl("bam_Oo|bam_2C|coverage_Oo|coverage_2C|knownGenes|repeatMasker|countFunction|filterRangesByFPKM|coverageForSummedPlot|positionCoverageDF|plotCummulativeFPKMBin", ls())])

################################################################################## functions
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
  
  # add rownames to new column 
  genomic_filtered_df %<>%
    mutate(feature_name = rownames(.)) %>% 
    set_rownames(NULL)
  
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
    select(-1) %>% 
    mutate(s_Oo = s_Oo / (number_of_reads["s_Oo"] * (tile_width / 1000)), 
           s_2C = s_2C / (number_of_reads["s_2C"] * (tile_width / 1000)))
  
  # output .csv
  all_coverage_counts_df_sum_wide %>% 
    mutate(pos = c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width))) %>% 
    write_csv(path = paste0("Park2013_", choosen_element_name, "_cummulativeFPKM_bins", tile_width, ".csv"), col_names = T)
  
  # long format for bin plot using tidyr::gather
  # syntax:
  # - gather all data except position (-position)
  # - name new key column stage (key = stage)
  # - name new value column fpkm (value = fpkm)
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

################################################################################## reading data
# setting working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/Park_2013_mm9")

# making TxDb object from knownGene gtf from UCSC, getting FPKM
knownGenes <- 
  makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm9/UCSC_knownGene_mm9_20170205.gtf.gz") %>% 
  genes(.)

# repeatMasker table
repeatMasker <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm9/UCSC_repeatMasker_mm9_20170205_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, gene_id = repName) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# combine knownGenes and repeatMasker tables for filtering out elements
genomic_filter <- c(knownGenes, repeatMasker)

################################################################################## reading .bam, getting coverage
# .bam files paths
filenames <- file.path(c("/common/DB/vfranke/Base/AccessoryData/Park_2013_GenDev/Merged/s_Oo.SE/s_Oo.SE.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Park_2013_GenDev/Merged/s_1C.SE/s_1C.SE.bam"))

# .bam files
bam_Oo <- readGAlignments(filenames[1])
bam_2C <- readGAlignments(filenames[2])

# coverage of BAM files
coverage_Oo <- coverage(bam_Oo)
coverage_2C <- coverage(bam_2C)

# getting library size in millions of reads
number_of_reads <- 
  read_lines("Park_2013_library_size.txt") %>% 
  as.integer() %>% 
  divide_by(10^6) %>% 
  set_names(c("s_Oo", "s_2C"))

################################################################################### reading full list, filtering by hand
# reading list of elements with expression in 2Cell
# MT2 and ORR1A0 solo LTRs 
MT2_ORR1A0_all <- 
  read_delim("MT2_ORR1A0_solo_LTRs_orderedByExpressionIn2cell_Park_mm9.txt", delim = "\t") %>% 
  dplyr::rename(element_name = repName)

# MT2 full length elements 
MT2_full <- 
  read_delim("MT2_full_orderedByExpressionIn2cell_Park_mm9.txt", delim = "\t") %>% 
  dplyr::select(-id) %>% 
  mutate(element_name = "MT2_full")

# combining all elements to one variable
all_elements <- 
  rbind(MT2_ORR1A0_all, MT2_full) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           element_name)) %>% 
  dplyr::filter(!(fullName %in% c("chr18:90631939-90632457:+|MT2_Mm"))) # filters out elements which go out of range of chromosome when expanded
  
################################################################################### finding overlaps between original repeat ranges and expanded ranges
# making original GRanges
all_elements_original_gr <- makeGRangesFromDataFrame(all_elements, keep.extra.columns = T)

# expanding ranges by 150kb in both directions, making expanded GRanges
all_elements_expanded_gr <- 
  all_elements %>% 
  mutate(end = start + 1.5e5, 
         start = start - 1.5e5) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# finding overlaps 
all_elements_overlaps <- findOverlaps(all_elements_expanded_gr, all_elements_original_gr)

# removing element which have other element in their 150kb expanded ranges
all_elements_overlaping_names <- 
  cbind(all_elements_original_gr[subjectHits(all_elements_overlaps)]$fullName, 
        all_elements_expanded_gr[queryHits(all_elements_overlaps)]$fullName) %>% 
  as.data.frame(.) %>% 
  set_colnames(c("fullName_original", "fullName_expanded")) %>% 
  dplyr::filter(fullName_original != fullName_expanded) %$% 
  fullName_expanded %>% 
  unique() %>% 
  as.character()

# filtering from original table  
all_elements %<>% 
  dplyr::filter(!(fullName %in% all_elements_overlaping_names))

################################################################################### get top 100 ordered by FPKM for each class
# top 100 MT2/ORR1A0
all_elements_top100 <- 
  all_elements %>%
  group_by(element_name) %>%
  arrange(desc(s_2C_fpkm)) %>%
  dplyr::slice(1:100) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-s_2C_fpkm) 

#####################################################################################
##################################################################################### loop through element types

# plots coverage plot and bined cummulative FPKM coverage plot
for(choosen_element_name in unique(all_elements_top100$element_name)){
  
  ##################################################################################### getting elements of one type
  # original ranges
  element_original <- 
    all_elements_top100 %>% 
    dplyr::filter(element_name == choosen_element_name) %>% 
    dplyr::select(-width)
  element_original_ranges <- makeGRangesFromDataFrame(element_original, keep.extra.columns = T)
  element_width <- max(width(element_original_ranges))
  
  # expand ranges 200kb up and downstream
  element_expanded <- 
    element_original %>% 
    mutate(end = start + 2e5, 
           start = start - 2e5) 
  element_expanded_ranges <- makeGRangesFromDataFrame(element_expanded, keep.extra.columns = T)
  
  ##################################################################################### filtering with knownGenes and RepeatMasker
  # removing original element ranges from filter (knownGenes and RepeatMasker)   
  genomic_filter_without_elements <- GenomicRanges::setdiff(genomic_filter, element_original_ranges, ignore.strand = T)
  seqlevels(genomic_filter_without_elements, force = T) <- seqlevels(bam_2C)
  
  # removing all regions overlaping with filtered repeatMasker/knownGene tables 
  element_genomic_filtered_df <- 
    lapply(element_expanded_ranges, function(x) GenomicRanges::setdiff(x, genomic_filter_without_elements, ignore.strand = T)) %>% 
    set_names(element_original$fullName) %>% 
    lapply(., as.data.frame) %>% 
    do.call(rbind, .) 
  element_genomic_filtered_ranges <- makeGRangesFromDataFrame(element_genomic_filtered_df)

  ##################################################################################### filtering elements above FPKM limit
  # getting counts of all samples with summarizeOverlaps function
  count_list <- lapply(list(bam_Oo, bam_2C), countFunction, input_features = element_genomic_filtered_ranges)
  
  # calculating fpkm from counts
  element_fpkm_df <- 
    do.call(cbind, count_list) %>% 
    set_colnames(c("s_Oo", "s_2C")) %>% 
    mutate(feature_names = rownames(.)) %>% 
    set_rownames(NULL) %>% 
    cbind(element_genomic_filtered_df, .) %>% 
    set_rownames(NULL) %>% 
    mutate(s_Oo = (s_Oo / (number_of_reads["s_Oo"] * (width / 1000))), 
           s_2C = (s_2C / (number_of_reads["s_2C"] * (width / 1000))))

  # removing element positions from FPKM data.frame 
  element_fpkm_ranges <- makeGRangesFromDataFrame(element_fpkm_df, keep.extra.columns = T)
  element_fpkm_ranges <- element_fpkm_ranges[-queryHits(findOverlaps(element_fpkm_ranges, element_original_ranges))]
  element_fpkm_filtered_df <- element_fpkm_df[element_fpkm_df$feature_names %in% element_fpkm_ranges$feature_names, ]
  
  # filtering ranges by FPKM values
  element_fpkm_filtered_grList <- 
    lapply(names(number_of_reads), 
           FUN = filterRangesByFPKM, 
           fpkm_limit = 1, 
           fpkm_filtered_df = element_fpkm_filtered_df, 
           genomic_filtered_df = element_genomic_filtered_df) %>% 
    set_names(names(number_of_reads))
  
  ##################################################################################### 
  # getting coverage 
  element_coverage_Oo <- lapply(X = 1:length(element_fpkm_filtered_grList$s_Oo), 
                                function(X) coverageForSummedPlot(X, 
                                                                  sample_coverage = coverage_Oo, 
                                                                  sample_name = "s_Oo", 
                                                                  element_grList = element_fpkm_filtered_grList)) 
  
  element_coverage_2C <- lapply(X = 1:length(element_fpkm_filtered_grList$s_2C), 
                                function(X) coverageForSummedPlot(X, 
                                                                  sample_coverage = coverage_2C, 
                                                                  sample_name = "s_2C", 
                                                                  element_grList = element_fpkm_filtered_grList)) 
  
  ##################################################################################### 
  # summing all stages coverage counts
  element_coverage_all_df <- 
    lapply(list(element_coverage_Oo, element_coverage_2C), positionCoverageDF) %>% 
    bind_rows(.) %>% 
    mutate(stage = rep(c("s_Oo", "s_2C"), each = (2 * 1.5e5 + element_width + 1)), 
           stage = factor(stage, levels = c("s_Oo", "s_2C")), 
           count = replace(count, count == 0, NA), 
           id = rep(1:(nrow(.)/2), length(number_of_reads)))
  
  # # plotting coverage
  # ggplot(element_coverage_all_df, aes(x = pos, y = count)) +
  #   geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
  #   scale_fill_manual(values = c("grey", "black")) +
  #   scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) + 
  #   coord_cartesian(ylim = c(0, 100)) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   facet_grid(stage ~ .) +
  #   ggsave(filename = paste0("Park2013_", choosen_element_name, "_summedCoverage.pdf"), width = 30, height = 10)
  # 
  ##################################################################################################### bin plot
  # plotCummulativeFPKMBin(merged_coverage_df = element_coverage_all_df, tile_width = 10000)
  
  
  ##################################################################################### Mann-Whitney test
  all_coverage_counts_df_merged <- element_coverage_all_df
  tile_width <- 10000
  
  # library size data.frame
  library_size_df <- 
    data.frame(library_size = as.integer(read_lines("Park_2013_library_size.txt"))) %>% 
    mutate(sample_name = c("s_Oo", "s_2C"), 
           library_size = library_size / 10^6) %>% 
    dplyr::select(2:1)
  
  # set tile width
  tile_width <- 10000
  
  # normalize to fpkm
  all_coverage_fpkm <- 
    all_coverage_counts_df_merged %>%
    dplyr::select(pos, fpkm = count, sample_name = stage) %>% 
    left_join(library_size_df, by = "sample_name") %>% 
    mutate(fpkm = fpkm / (library_size * (tile_width / 1000)), 
           sample_name = factor(sample_name, levels = library_size_df$sample_name)) %>% 
    dplyr::select(sample_name, pos, fpkm) %>% 
    tidyr::spread(key = sample_name, value = fpkm)
  
  # divide to upstream/downstream bins
  all_coverage_fpkm_binned <- 
    rbind(all_coverage_fpkm %>% #upstream
            filter(pos >= -150000 & pos < 0) %>%
            mutate(bin = rev(gl(ceiling(n() / tile_width), tile_width, n())), 
                   stream = "upstream", 
                   pos = abs(pos)),
          all_coverage_fpkm %>% #downstream
            filter(pos > element_width & pos <= (150000 + element_width)) %>%
            mutate(bin = gl(ceiling(n() / tile_width), tile_width, n()), 
                   stream = "downstream", 
                   pos = pos - element_width)) %>% 
    mutate(bin = as.integer(bin))
  
  ######## compare 2C upstream/downstream
  # take individual stage data and do 
  wilcox_results <- 
    all_coverage_fpkm_binned %>% 
    dplyr::select(pos, bin, stream, fpkm = s_2C) %>% 
    reshape2::dcast(pos + bin ~ stream, value.var = "fpkm") %>% 
    group_by(bin) %>% 
    do(tidy(wilcox.test(.$downstream, .$upstream, alternative = "greater", paired = F))) %>% 
    ungroup(wilcox_results) %>% 
    as.data.frame() %>% 
    mutate(significant = p.value <= 0.05) %>% 
    dplyr::select(-method)
  
  ######## compare downstream oocyte/2C
  wilcox_results <- 
    all_coverage_fpkm_binned %>% 
    dplyr::select(pos, bin, stream, s_2C, s_Oo) %>% 
    dplyr::filter(stream == "downstream") %>% 
    dplyr::select(-stream) %>% 
    group_by(bin) %>% 
    do(tidy(wilcox.test(.$s_2C, .$s_Oo, alternative = "greater", paired = F))) %>% 
    ungroup(wilcox_results) %>% 
    as.data.frame() %>% 
    mutate(significant = p.value <= 0.05) %>% 
    dplyr::select(-method)
  
}

