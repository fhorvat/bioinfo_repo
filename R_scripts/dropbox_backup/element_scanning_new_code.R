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

# # removes variables from enviorment (all except loaded .bam files, coverages and functions)
# loaded_functions <- intersect(ls(), lsf.str())
# loaded_variables <- grep("s_.*_bam|s_.*_coverage", ls(), value = T)
# rm(list = ls()[!(ls() %in% c(loaded_functions, loaded_variables)) | grepl("element", ls())])

################################################################################## functions
countFunction <- function(feature, bam_reads, single_end){
  
  register(MulticoreParam(workers = 8))
  se_filtered <- summarizeOverlaps(features = feature,
                                   reads = bam_reads, 
                                   mode = "Union", 
                                   singleEnd = single_end, 
                                   ignore.strand = TRUE)
  
  return(as.data.frame(assay(se_filtered)))
  
}

filterRangesByFPKM <- function(sample_name, element_fpkm, element_ranges, fpkm_limit){
  
  # get names of features above FPKM limit
  fpkm_df_names_limit <- element_fpkm[element_fpkm[, sample_name] > fpkm_limit, "feature_name"]
  
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
  
  # reduce sample coverage to chromosome where sample is, create RangesList for sample and also reduce it to one chromosme  
  sample_coverage_reduced <- sample_coverage[seqnames(element_gr[1])]
  element_rangesList <- as(element_gr, "RangesList")[seqnames(element_gr[1])]
  
  # find coverage
  element_coverage <- Views(sample_coverage_reduced, element_rangesList)
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
  
  # getting all unique position across all elements in one sample
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
  # also normalize for library size = RPM/FPM
  all_coverage_counts_df_sum_wide <- 
    merged_coverage_df %>% 
    left_join(library_size_df, by = "sample") %>% 
    mutate(count = count / library_size) %>% 
    dplyr::select(sample, pos, count) %>% 
    mutate(count = replace(count, is.na(count), 0)) %>% 
    tidyr::spread(key = sample, value = count) 
  
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
            dplyr::select(indx, starts_with("s_"), position),
          all_coverage_counts_df_sum_wide %>% 
            dplyr::filter(pos > element_width & pos <= (150000 + element_width)) %>% 
            group_by(indx = gl(ceiling(nrow(.)/tile_width), tile_width, nrow(.))) %>%
            dplyr::select(-pos) %>%
            summarise_each(funs(sum)) %>% 
            mutate(position = paste(seq(0, (150000 - tile_width), tile_width), "downstream", sep = "_"))) %>%
    select(-1) %>% 
    mutate_at(.cols = vars(starts_with("s_")), 
              .funs = funs(divide_by(., (tile_width / 1000))))
  
  # output .csv
  all_coverage_counts_df_sum_wide %>% 
    mutate(pos = c(seq(-150000, -tile_width, tile_width), 
                   0, 
                   seq(element_width, ((150000 + element_width) - tile_width), tile_width))) %>% 
    write_csv(path = paste0("./results/", file_name, "_cummulativeFPKM_", tile_width, "_bins.csv"), col_names = T)
  
  # get sample order based on summed FPKM 
  sample_order <- 
    all_coverage_counts_df_sum_wide %>% 
    summarise_at(.cols = vars(starts_with("s_")), sum) %>% 
    tidyr::gather(key = sample, value = sum_fpkm) %>% 
    arrange(desc(sum_fpkm)) 
  
  # long format for bin plot using tidyr::gather
  # syntax:
  # - gather all data except position (-position)
  # - name new key column sample (key = sample)
  # - name new value column fpkm (value = fpkm)
  all_coverage_counts_df_sum_long <- 
    all_coverage_counts_df_sum_wide %>% 
    tidyr::gather(key = sample, value = fpkm, -position) %>% 
    mutate(element = ifelse(position == "element", "in_element",  "out_element"), 
           width = ifelse(position == "element", element_width, tile_width), 
           pos = rep(c(seq(-150000, -tile_width, tile_width), 0, seq(element_width, ((150000 + element_width) - tile_width), tile_width)), nrow(library_size_df)), 
           sample = factor(sample, levels = sample_order$sample)) %>% 
    dplyr::select(pos, fpkm, element, sample, position, width) %>% 
    arrange(sample) %>% 
    mutate(sample = factor(sample, levels = library_size_df$sample))
  
  ################################################################################## plot
  # get shades of grey for color pallete, change the color of highest value to orange
  cols <- 
    colorRampPalette(brewer.pal(n = nrow(library_size_df), name = "Greys"))(nrow(library_size_df)) %>% 
    set_names(sample_order$sample) 
  cols[sample_order$sample[1]] <- "orange"
  
  # bin plot
  ggplot(all_coverage_counts_df_sum_long, aes(x = pos, y = fpkm)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + width, ymax = fpkm, fill = sample)) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width), 
                       breaks = c(seq(-150000, -10000, 10000), 0, seq(element_width, (140000 + element_width), 10000)), 
                       labels = c(seq(150, 10, -10), element_name, seq(0, 140, 10)), 
                       name = "bin") + 
    scale_y_continuous(name = "Cummulative FPKM") +  
    coord_cartesian(ylim = c(0, 2000)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = paste0("./results/", file_name, "_cummulativeFPKM_", tile_width, "_binPlot.pdf"), width = 30, height = 10)
  
}

################################################################################## reading data
# setting working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/cow_expression/Fugaku_test")

# .bam files paths
filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))

# .bam files
s_GV_bam <- readGAlignmentPairs(filenames[1])
s_1C_bam <- readGAlignmentPairs(filenames[2])
s_2C_bam <- readGAlignmentPairs(filenames[3])
s_2Ca_bam <- readGAlignmentPairs(filenames[4])
s_4C_bam <- readGAlignmentPairs(filenames[5])

# coverage of BAM files
s_GV_coverage <- coverage(s_GV_bam)
s_1C_coverage <- coverage(s_1C_bam)
s_2C_coverage <- coverage(s_2C_bam)
s_2Ca_coverage <- coverage(s_2Ca_bam)
s_4C_coverage <- coverage(s_4C_bam)

# get library size in million of reads
library_size_df <- 
  sapply(X = 1:length(logs_filenames), 
         FUN = function(X){
           read_delim(logs_filenames[X], delim = "\t", col_names = F)[8, 2] %>%
             as.integer(.)
         }) %>% 
  data.frame(library_size = . / 10^6) %>% 
  mutate(sample = factor(c("s_GV", "s_1C", "s_2C", "s_2Ca", "s_4C"), levels = c("s_GV", "s_1C", "s_2C", "s_2Ca", "s_4C"))) %>% 
  dplyr::select(sample, library_size) 

################################################################################## reading annotations 
# making TxDb object from refseq gtf
genes_ranges <- 
  makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz") %>% 
  genes(.)

# repeatMasker table
repeatMasker <- 
  read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, repClass, repFamily) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T) 

######################################################################### getting elements
# read from table, remove 
all_elements_original <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_full_orderedByFPKMIn2cell.txt", delim = "\t") %>% 
  mutate(repName = "MT2_full", 
         fullName = paste0(seqnames, ":", 
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           repName)) %>% 
  dplyr::select(-c(id, width, repName)) %>% 
  dplyr::select(1:4, 6:5)
  
# make GRanges
all_elements_original_ranges <- makeGRangesFromDataFrame(all_elements_original, keep.extra.columns = T)

##################################################################################### filtering based on distance to nearest gene  
# find distance to nearest gene 
all_elements_original_ranges$distance_to_genes <- mcols(distanceToNearest(all_elements_original_ranges, genes_ranges, ignore.strand = T))$distance

experiment <- "Fugaku"
element_name <- "MT2_full"
distances <- c(0, 50000)

##################################################################################### loop over distance to nearest gene

for(min_distance in distances){
  
  # filter original ranges
  element_original_ranges <- all_elements_original_ranges[all_elements_original_ranges$distance_to_genes >= min_distance]
  element_original_ranges <- element_original_ranges[order(element_original_ranges$FPKM, decreasing = T)]
  
  # take top 100, remove RPKM data, get max.element width
  if(length(element_original_ranges) >= 100){
    element_original_ranges <- element_original_ranges[1:100]
  }
  
  mcols(element_original_ranges)$FPKM <- NULL
  element_width <- max(width(element_original_ranges))
  
  # expand ranges 200kb up and downstream
  element_expanded_ranges <- 
    as.data.frame(element_original_ranges) %>% 
    mutate(end = start + 2e5, 
           start = start - 2e5) %>% 
    makeGRangesFromDataFrame(., keep.extra.columns = T)
  
  ##################################################################################### filtering with repeatMasker and/or genes
  # removing original element ranges from filter
  genomic_filter <- c(repeatMasker)
  genomic_filter <- GenomicRanges::setdiff(genomic_filter, element_original_ranges, ignore.strand = T)
  seqlevels(genomic_filter, force = T) <- seqlevels(s_2C_bam)
  
  # setting file name
  file_name <- paste0(experiment, "_", element_name, "_dist", (min_distance / 1000), "kb_top", length(element_original_ranges))
  
  # filtering with genomic filter, getting data.frame
  element_genomic_filtered_df <- 
    lapply(element_expanded_ranges, function(x) GenomicRanges::setdiff(x, genomic_filter, ignore.strand = T)) %>% 
    set_names(element_expanded_ranges$fullName) %>% 
    lapply(., as.data.frame) %>% 
    do.call(rbind, .) %>% 
    mutate(feature_name = rownames(.)) %>% 
    set_rownames(NULL)
  
  ##################################################################################### filtering with FPKM threshold 
  # getting counts of all samples with summarizeOverlaps, joining with library size, calculating FPKM
  element_genomic_filtered_fpkm <- 
    lapply(list(s_GV_bam, s_1C_bam, s_2C_bam, s_2Ca_bam, s_4C_bam), 
           FUN = countFunction, 
           feature = makeGRangesFromDataFrame(element_genomic_filtered_df),
           single_end = FALSE) %>% 
    bind_cols(.) %>% 
    set_colnames(library_size_df$sample) %>% 
    cbind(element_genomic_filtered_df, .) %>% 
    dplyr::select(-strand) %>% 
    tidyr::gather(key = sample, value = fpkm, -c(1:5)) %>% 
    mutate(sample = factor(sample, levels = library_size_df$sample)) %>% 
    left_join(library_size_df, by = "sample") %>% 
    mutate(fpkm = fpkm / (library_size * (width / 1000)), 
           sample = factor(sample, levels = library_size_df$sample)) %>% 
    dplyr::select(-library_size) %>% 
    tidyr::spread(key = sample, value = fpkm)
  
  # removing MT2 positions from FPKM data.frame 
  element_genomic_filtered_fpkm_gr <- makeGRangesFromDataFrame(element_genomic_filtered_fpkm, keep.extra.columns = T)
  element_genomic_filtered_fpkm_gr <- element_genomic_filtered_fpkm_gr[-queryHits(findOverlaps(element_genomic_filtered_fpkm_gr, element_original_ranges))]
  element_genomic_filtered_fpkm <- element_genomic_filtered_fpkm[element_genomic_filtered_fpkm$feature_name %in% element_genomic_filtered_fpkm_gr$feature_name, ]
  
  # filtering ranges by FPKM values 
  # - returns GRList with FPKM filtered element ranges for each .bam file
  element_fpkm_filtered <-
    lapply(as.character(library_size_df$sample), 
           FUN = filterRangesByFPKM, 
           element_fpkm = element_genomic_filtered_fpkm, 
           element_ranges = element_genomic_filtered_df,
           fpkm_limit = 1) %>% 
    set_names(library_size_df$sample)
  
  #########################################################################################################
  # getting coverage as list of coverage vectors (with names = relative position) for each of 100 elements (6x, for each sample/.bam file)
  # GV
  s_GV_element_coverage <- lapply(X = 1:length(element_fpkm_filtered[["s_GV"]]), 
                                  function(X) coverageForSummedPlot(X, 
                                                                    sample_coverage = s_GV_coverage, 
                                                                    sample_name = "s_GV", 
                                                                    element_grList = element_fpkm_filtered)) 
  
  # 1-cell
  s_1C_element_coverage <- lapply(X = 1:length(element_fpkm_filtered[["s_1C"]]), 
                                   function(X) coverageForSummedPlot(X, 
                                                                     sample_coverage = s_1C_coverage, 
                                                                     sample_name = "s_1C", 
                                                                     element_grList = element_fpkm_filtered)) 
  
  # 2-cell
  s_2C_element_coverage <- lapply(X = 1:length(element_fpkm_filtered[["s_2C"]]), 
                                  function(X) coverageForSummedPlot(X, 
                                                                    sample_coverage = s_2C_coverage, 
                                                                    sample_name = "s_2C", 
                                                                    element_grList = element_fpkm_filtered)) 
  
  # 2-cell + aphicolidin 
  s_2Ca_element_coverage <- lapply(X = 1:length(element_fpkm_filtered[["s_2Ca"]]), 
                                  function(X) coverageForSummedPlot(X, 
                                                                    sample_coverage = s_2Ca_coverage, 
                                                                    sample_name = "s_2Ca", 
                                                                    element_grList = element_fpkm_filtered)) 
  
  # 4-cell
  s_4C_element_coverage <- lapply(X = 1:length(element_fpkm_filtered[["s_4C"]]), 
                                   function(X) coverageForSummedPlot(X, 
                                                                     sample_coverage = s_4C_coverage, 
                                                                     sample_name = "s_4C", 
                                                                     element_grList = element_fpkm_filtered)) 
  
  ######################################################################################################### 
  # summing coverage of all elements of one sample
  element_coverage_all_df <- 
    lapply(list(s_GV_element_coverage, s_1C_element_coverage, s_2C_element_coverage, s_2Ca_element_coverage, s_4C_element_coverage), positionCoverageDF) %>% 
    bind_rows(.) %>% 
    mutate(sample = rep(library_size_df$sample, each = (2 * 1.5e5 + element_width + 1)), 
           sample = factor(sample, levels = library_size_df$sample), 
           count = replace(count, count == 0, NA))
  # id = rep(1:(nrow(.)/2), nrow(library_size_df))
  
  # plotting coverage
  ggplot(element_coverage_all_df, aes(x = pos, y = count)) +
    geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
    scale_fill_manual(values = c("grey", "black")) +
    scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) + 
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_grid(sample ~ .) +
    ggsave(filename = paste0("./results/", file_name, "_summedCoverage.pdf"), width = 30, height = 10)
  
  ##################################################################################################### bin plot
  plotCummulativeFPKMBin(merged_coverage_df = element_coverage_all_df, tile_width = 10000)
  
  cat(file_name, "done", "\n")
  
}  
