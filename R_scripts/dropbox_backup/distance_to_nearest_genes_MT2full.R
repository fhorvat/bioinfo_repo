library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(cowplot)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(BiocParallel)

options(bitmapType = "cairo")

################################################################################## reading data
# setting working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/distance_to_nearest_genes")

# reading original element lists (MT2 and ORR1A0 solo LTRs ordered by expression in 2Cell)
MT2_ORR1A0_all <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_ORR1A0_solo_LTRs_orderedByFPKMin2cell.txt", delim = "\t") %>% 
  dplyr::select(-4)

MT2_full <- 
  read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_full_orderedByFPKMIn2cell.txt", delim = "\t") %>% 
  dplyr::select(-c(4, 8)) %>% 
  mutate(repName = "MT2_full")

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
                                          "chr13:76077210-76083662:-|MT2_full")

all_elements_filtered <- subset(all_elements_filtered, !(all_elements_filtered$fullName %in% c(hand_picked_filter_elements_MT2_solo, 
                                                                                               hand_picked_filter_elements_ORR1A0, 
                                                                                               hand_picked_filter_elements_MT2_full)))
# top 100 MT2/ORR1A0
all_elements_original <- 
  all_elements_filtered %>%
  group_by(repName) %>%
  arrange(desc(FPKM)) %>%
  dplyr::slice(1:100) %>% 
  ungroup() %>% 
  dplyr::select(-FPKM) %>% 
  as.data.frame(.)

# making TxDb object from knownGene gtf from UCSC, getting FPKM
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz")
knownGenes_gtf_gr <- genes(knownGenes_gtf)

##################################################################################### run from here
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

##################################################################################### 
# find distance between MT2 full elements and knownGenes
distance_to_nearest <- distanceToNearest(MT2OriginalGR, knownGenes_gtf_gr, ignore.strand = T)

# extract nearest genes
distance_to_nearest_genes <- knownGenes_gtf_gr[subjectHits(distance_to_nearest)]

##################################################################################### 
# extract nearest MT2s - 
distance_to_nearest_elements <- 
  MT2OriginalGR[queryHits(distance_to_nearest)] %>% 
  as.data.frame(.) %>% 
  mutate(distance = mcols(distance_to_nearest)$"distance", 
         nearest_gene_start = start(distance_to_nearest_genes), 
         gene_id = mcols(distance_to_nearest_genes)$gene_id) %>% 
  dplyr::filter(distance <= 150000) %>% 
  mutate(relative_distance = start - nearest_gene_start, 
         gene_stream = ifelse(relative_distance < 0, "down", "up")) %>% 
  dplyr::select(-c(relative_distance, nearest_gene_start)) %>% 
  mutate(distance = ifelse(gene_stream == "up", -distance, distance))

# adding ones overlaping LTRs twice (for both upstream and downstream class)
distance_to_nearest_elements_duplicated <- 
  distance_to_nearest_elements %>% 
  dplyr::filter(distance == 0) %>% 
  mutate(gene_stream = "down")

distance_to_nearest_elements <- 
  bind_rows(distance_to_nearest_elements, 
            distance_to_nearest_elements_duplicated) 

distance_to_nearest_elements %<>% 
  dplyr::filter(distance != 0)

median(abs(distance_to_nearest_elements$distance))
nrow(distance_to_nearest_elements)
##################################################################################### 
# plot histogram
plot_hist <- 
  ggplot(data = distance_to_nearest_elements, aes(x = distance)) +
  geom_histogram(binwidth = 1000)

# boxplot downstream
plot_box <- 
  ggplot(data = distance_to_nearest_elements, aes(x = gene_stream, y = distance)) +
  geom_boxplot(fill = NA) +
  geom_jitter(color = "red", size = 1, height = 0.05) +
  coord_flip()

# plot in grid
plot_grid(plot_hist, plot_box, nrow = 2, align = "v") +
  ggsave("bins_1000_noOverlap_2cExpressed.pdf", width = 20, height = 10)

