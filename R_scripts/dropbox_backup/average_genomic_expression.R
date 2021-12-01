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

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/statistics")

################################################################################## functions
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
  countsZero <- rep(0, 4e5 + 1)
  names(countsZero) <- MT2SampleExpanded$start : MT2SampleExpanded$end
  MT2SampleCoverageCountsFull <- replace(countsZero, names(countsZero) %in% names(MT2SampleCoverageCounts), MT2SampleCoverageCounts)
  
  MT2SampleCoveragePosition <- as.numeric(names(MT2SampleCoverageCountsFull))
  MT2SampleCoveragePosition <- MT2SampleCoveragePosition - MT2SampleExpanded$start - 2e5
  
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
  position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000"):which(rownames(position_matrix) == as.character(150000 + element_width)), ]
  
  allMT2CoverageCountsSummedFinal <- rowSums(position_matrix)
  allMT2CoverageCountsSummedFinal <- data.frame(pos = as.numeric(names(allMT2CoverageCountsSummedFinal)), 
                                                count = allMT2CoverageCountsSummedFinal)
  allMT2CoverageCountsSummedFinal$element <- ifelse((allMT2CoverageCountsSummedFinal$pos >= 0 & allMT2CoverageCountsSummedFinal$pos <= element_width),
                                                    "in_element", 
                                                    "out_element")
  allMT2CoverageCountsSummedFinal$stage <- stage
  return(allMT2CoverageCountsSummedFinal)
}

# removes outliers from vector
removeOutliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

################################################################################## reading annotations 
# get library size in million of reads
library_size_df <- 
  read_delim(file = "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/cow_expression/Graf_2014_library_size.txt", delim = "\t") %>% 
  mutate(sample = factor(sample, levels = sample))

# making TxDb object from genes gtf, getting genes 
genes_ranges_list <- 
  makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz") %>% 
  genes(.) %>% 
  GenomicRanges::split(., seqnames(.))

# # repeatMasker table
# repeatMasker <- 
#   read_delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", delim = "\t") %>% 
#   dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName) %>% 
#   makeGRangesFromDataFrame(., keep.extra.columns = T) 

# get 100kb binned chromosomes
binned_chr_list <- 
  rtracklayer::import.bed("/common/WORK/fhorvat/reference/mouse/mm10/mm10_genome_binned_100kb.bed") %>% 
  GenomicRanges::split(., seqnames(.))

################################################################################## mask knownGenes from binned ranges
# filter out smaller scaffolds 
genes_ranges_list <- genes_ranges_list[!str_detect(string = names(genes_ranges_list), pattern = "random|Un")] 
binned_chr_list <- binned_chr_list[!str_detect(string = names(binned_chr_list), pattern = "random|Un")]

# order 
genes_ranges_list <- genes_ranges_list[order(names(genes_ranges_list))]
binned_chr_list <- binned_chr_list[order(names(binned_chr_list))]

# filter knownGenes from binned chromosomes
binned_chr_list_masked <- lapply(X = 1:length(binned_chr_list), FUN = function(X){
  
  # get all genes on one chromosome
  one_chr_genes <- genes_ranges_list[[X]]
  
  # get all bins on one chromosome
  one_chr_binned <- binned_chr_list[[X]] 
  
  # remove ranges belonging to genes from binned chromosome ranges 
  one_chr_binned_masked <- 
    lapply(one_chr_binned, function(x) GenomicRanges::setdiff(x, one_chr_genes, ignore.strand = T)) %>% 
    set_names(., paste0(names(binned_chr_filtered)[X], "_bin", 1:length(.))) %>% 
    GRangesList(.)

})

# merge GRangesLists 
binned_chr_list_masked_merged <- do.call(c, binned_chr_list_masked)

################################################################################## 
# .bam files paths
filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAm.bam",
                         "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WE.bam"))

# reading .bam files
bam_2C <- readGAlignmentPairs(filenames[3])

# finding coverage of BAM files
coverage_2C <- coverage(bam_2C)

################################################################################## 
# get expression in 2C over filtered binned ranges
# bamfiles <- BamFileList("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
#                         yieldSize = 2000000)
# se <- summarizeOverlaps(features = binned_chr_list_masked_merged,
#                         reads = bamfiles, 
#                         mode = "Union", 
#                         singleEnd = FALSE, 
#                         ignore.strand = TRUE)
# 
# # get counts in 2C
# binned_ranges_counts <- 
#   as.data.frame(assay(se)) %>% 
#   set_colnames("counts") %>% 
#   mutate(bin_name = rownames(.)) %>% 
#   set_rownames(NULL) %>% 
#   cbind(as.data.frame(unlist(binned_chr_list), row.names = NULL), .) %T>% 
#   write_csv("counts_2C_on_100kb_binned_genome.csv")

# get counts on 100kb binned chromosomes with masked genes
binned_ranges_counts <- read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/statistics/counts_2C_on_100kb_binned_genome.csv")
binned_ranges_counts_ordered <- 
  binned_ranges_counts %>% 
  arrange(desc(counts))

# get mean count
mean_genomic_count <- 
  binned_ranges_counts$counts %>% 
  removeOutliers %>% 
  na.omit %>% 
  mean

#################################################################################
# making TxDb object from knownGene gtf from UCSC, getting FPKM
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz")
knownGenes_gtf_gr <- genes(knownGenes_gtf)

# reading repeatMasker table, filtering for MT/MT2/ORR
rptmsk <- read.delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMasker_mm10_20170209_all_fields.txt.gz", stringsAsFactors = F, header = T) 
rptmsk <- rptmsk[, c("genoName", "genoStart", "genoEnd", "strand", "repName")]
colnames(rptmsk) <- c("seqnames", "start", "end", "strand", "gene_id")
rptmsk <- makeGRangesFromDataFrame(rptmsk, keep.extra.columns = T)

rptmskFiltered <- rptmsk[grep("MT|ORR", mcols(rptmsk)$gene_id)]
rptmskFiltered <- rptmskFiltered[-grep("MMTV-int", mcols(rptmskFiltered)$gene_id)]

knownGenesAndRepeatMasker <- c(knownGenes_gtf_gr, rptmsk)
all_elements_original <- as.data.frame(all_elements_original)

#################################################################################
#################################################################################
# read top100 MT2 full elements used in original paper figure
MT2_top <- read_csv("MT2full_originalRanges_top100in2cell_figure6BinOriginalPaper.csv")

# make original element ranges
MT2OriginalGR <- makeGRangesFromDataFrame(MT2_top, keep.extra.columns = T)

# expand ranges 200kb up and downstream
MT2Expanded <- MT2_top
MT2Expanded$end <- MT2Expanded$start + 2e5 
MT2Expanded$start <- MT2Expanded$start - 2e5 
MT2ExpandedGR <- makeGRangesFromDataFrame(MT2Expanded, keep.extra.columns = T)

# joining filtered knownGenes and repeatMasker tables (to be used later for filtering reads)
knownGenesAndRepeatMaskerFiltered <- GenomicRanges::setdiff(knownGenesAndRepeatMasker, MT2OriginalGR, ignore.strand = T)
seqlevels(knownGenesAndRepeatMaskerFiltered, force = T) <- seqlevels(bam_2C)

# filtering with knownGenes and RepeatMasker
MT2ExpandedFilteredGRList <- lapply(MT2ExpandedGR, function(x) GenomicRanges::setdiff(x, knownGenesAndRepeatMaskerFiltered, ignore.strand = T))
names(MT2ExpandedFilteredGRList) <- MT2_top$repNames 
MT2ExpandedFilteredGRList <- MT2ExpandedFilteredGRList[!names(MT2ExpandedFilteredGRList) %in% MT2_hand_picked_filter]

#################################################################################
# getting FPKM values for filtered features
MT2_filtered_ranges <- lapply(MT2ExpandedFilteredGRList, as.data.frame)
MT2_filtered_ranges_df <- do.call(rbind, MT2_filtered_ranges)
MT2_filtered_ranges_gr <- makeGRangesFromDataFrame(MT2_filtered_ranges_df)

### calculating fpkm values
# getting library size in millions of reads
logs_filenames <- file.path("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out") 

number_of_reads <- sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- "2C"
number_of_reads <- number_of_reads / 10^6

# getting counts of all samples with summarizeOverlaps function
count_list <- lapply(list(bam_2C), count_function)
count_df <- do.call(cbind, count_list)
colnames(count_df) <- "bam_2C"
count_df <- as.data.frame(count_df)
count_df$feature_names <- rownames(count_df)
rownames(count_df) <- NULL
MT2_filtered_ranges_df_counts <- cbind(MT2_filtered_ranges_df, count_df)
rownames(MT2_filtered_ranges_df_counts) <- NULL

# calculating fpkm from counts
MT2_filtered_ranges_df_fpkm <- MT2_filtered_ranges_df_counts
MT2_filtered_ranges_df_fpkm[, "bam_2C"] <- MT2_filtered_ranges_df_fpkm$bam_2C / 
  (number_of_reads["2C"] * (MT2_filtered_ranges_df_fpkm$width / 1000))

# removing MT2 positions from fpkm 
MT2_filtered_ranges_df_fpkm_gr <- makeGRangesFromDataFrame(MT2_filtered_ranges_df_fpkm, keep.extra.columns = T)
MT2_filtered_ranges_df_fpkm_gr_overlaps <- findOverlaps(MT2_filtered_ranges_df_fpkm_gr, MT2OriginalGR)
MT2_filtered_ranges_df_fpkm_gr <- MT2_filtered_ranges_df_fpkm_gr[-queryHits(MT2_filtered_ranges_df_fpkm_gr_overlaps)]
MT2_filtered_ranges_df_fpkm_filtered <- MT2_filtered_ranges_df_fpkm[MT2_filtered_ranges_df_fpkm$feature_names %in% MT2_filtered_ranges_df_fpkm_gr$feature_names, ]

#########################################################################################################
# filtering ranges by FPKM values
element_width <- max(width(MT2OriginalGR))

all_stages_bam_names <- "bam_2C"
MT2_filtered_ranges_fpkm_filtered_grList_all <- lapply(all_stages_bam_names, filterRangesByFPKM, fpkm_limit = 10000)
names(MT2_filtered_ranges_fpkm_filtered_grList_all) <- all_stages_bam_names

# getting coverage counts
allMT2CoverageCounts_2C <- lapply(X = 1:length(MT2_filtered_ranges_fpkm_filtered_grList_all$"bam_2C"), 
                                  function(X) coverageForSummedPlot(X, coverage_2C, "bam_2C")) 

#########################################################################################################
# summing all stages coverage counts
all_stages_counts <- list(allMT2CoverageCounts_2C)
all_stages_names <- "2C"
all_coverage_counts <- lapply(X = 1:length(all_stages_counts), 
                              function(X) positionCoverageDF(allMT2CoverageCounts = all_stages_counts[[X]], stage = all_stages_names[X]))
all_coverage_counts_df <- do.call(rbind, all_coverage_counts)
all_coverage_counts_df$stage <- factor(all_coverage_counts_df$stage, "2C")
all_coverage_counts_df$count[all_coverage_counts_df$count == 0] <- NA
all_coverage_counts_df$pos[is.na(all_coverage_counts_df$count)] <- NA

# plotting with ggplot scaled
ggplot(all_coverage_counts_df, aes(x = pos, y = count)) +
  geom_rect(aes(xmin = pos, ymin = 0, xmax = pos + 1, ymax = count, fill = element)) +
  scale_fill_manual(values = c("grey", "black")) +
  scale_x_continuous(limits = c(-1.5e5, 1.5e5 + element_width)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # facet_grid(stage ~ .) +
  ggsave(filename = paste0(element_name, "_noFPKM_filter_summedCoverage_.pdf"), width = 30, height = 10)


#########################################################################################################
# summing all stages coverage counts
allMT2CoverageCounts <- allMT2CoverageCounts_2C
stage <- "2C"

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
position_matrix <- position_matrix[which(rownames(position_matrix) == "-150000"):which(rownames(position_matrix) == as.character(150000 + element_width)), ]

# get data.frame from position matrix, filter out element, add info about stream
position_matrix_df <- 
  data.frame(position_matrix) %>% 
  mutate(pos = rownames(.)) %>% 
  set_rownames(NULL) %>% 
  dplyr::mutate(pos = as.numeric(pos)) %>% 
  dplyr::filter(!(pos >= 0 & pos <= element_width)) %>% 
  mutate(stream = ifelse(pos < 0, "upstream", "downstream"))

# take only positions downstream from element, add bins
position_matrix_df_downstream <- 
  position_matrix_df %>% 
  dplyr::filter(stream == "downstream") %>%
  dplyr::select(-stream) %>% 
  mutate(pos = pos - element_width, 
         bin = gl(ceiling(nrow(.) / 10000), 10000, nrow(.))) %>% 
  group_by(bin) %>% 
  summarise_at(.cols = vars(starts_with("X")), .funs = sum) %>% 
  tidyr::gather(key = element_number, value = count, -bin)

not_expressed_elements <- 
  position_matrix_df_downstream %>% 
  group_by(element_number) %>% 
  summarise(element_summed = sum(count)) %>% 
  dplyr::filter(element_summed == 0) %$%
  element_number
  
# get position of bin which is smaller than mean genomic expression
position_matrix_df_slice <- 
  position_matrix_df_downstream %>%
  dplyr::filter(!(element_number %in% not_expressed_elements) & (count < mean_genomic_count)) %>% 
  group_by(element_number) %>% 
  slice(1) %>% 
  as.data.frame()

