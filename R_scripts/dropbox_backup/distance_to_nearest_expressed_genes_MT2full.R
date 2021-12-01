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

# loaded_variables <- grep("s_.*_bam|s_.*_coverage", ls(), value = T)
# rm(list = ls()[!(ls() %in% loaded_variables) | grepl("element", ls())])

################################################################################## reading data
# setting working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/distance_to_nearest_genes/downstream_genes")

# making TxDb object from knownGene gtf from UCSC, getting FPKM of genes
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz")
knownGenes_gtf_gr <- genes(knownGenes_gtf)

# library size data.frame
logs_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out", 
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAmLog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WELog.final.out"))

library_size_df <- 
  data.frame(sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))) %>% 
  set_colnames("library_size") %>% 
  mutate(sample_name = c("s_GV", "s_1C", "s_2C", "s_2Ca", "s_4C"), 
         library_size = library_size / 10^6) %>% 
  dplyr::select(2:1)

################################################################################## expression of genes
# get expression of genes in 2cell
bamfiles <- BamFileList("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WE.bam", 
                        yieldSize = 2000000)
register(MulticoreParam(workers = 8))
se <- summarizeOverlaps(features = knownGenes_gtf_gr,
                        reads = bamfiles,
                        mode = "Union",
                        singleEnd = FALSE,
                        ignore.strand = TRUE)

# calculate fpkm 
knownGenes_expressed_id <- 
  as.data.frame(assay(se)) %>% 
  set_colnames("FPKM") %>% 
  mutate(gene_id = rownames(.), 
         library_size = library_size_df %>% filter(sample_name == "s_2C") %$% library_size, 
         gene_width = width(knownGenes_gtf_gr), 
         FPKM = (FPKM / (library_size * (gene_width / 1000)))) %>% 
  dplyr::select(gene_id, FPKM) %>% 
  dplyr::filter(FPKM > 0) %$%
  gene_id

# take only expressed genes
knownGenes_expressed <- knownGenes_gtf_gr[mcols(knownGenes_gtf_gr)$gene_id %in% knownGenes_expressed_id]

##################################################################################### run from here
# top 100 MT2 full elements used in original paper
MT2_top <- read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/statistics/MT2full_originalRanges_top100in2cell_figure6BinOriginalPaper.csv")
MT2OriginalGR <- makeGRangesFromDataFrame(MT2_top, keep.extra.columns = T)

##################################################################################### 
# finds distance of elements to nearest gene (all or expressed) + plot
distanceStatistics <- function(gene_type, plot){

  if(gene_type == "all"){
    genes_ranges <- knownGenes_gtf_gr
  }else{
    if(gene_type == "expressed"){
      genes_ranges <- knownGenes_expressed
    }
  }
  
  # find distance 
  distance_to_nearest <- distanceToNearest(MT2OriginalGR, genes_ranges, ignore.strand = T)
  distance_to_nearest_genes <- genes_ranges[subjectHits(distance_to_nearest)]
  
  # extract nearest MT2s - 
  distance_to_nearest_elements <- 
    MT2OriginalGR[queryHits(distance_to_nearest)] %>% 
    as.data.frame(.) %>% 
    mutate(distance = mcols(distance_to_nearest)$"distance", 
           nearest_gene_start = start(distance_to_nearest_genes), 
           gene_id = mcols(distance_to_nearest_genes)$gene_id) %>% 
    # dplyr::filter(distance <= 150000) %>%
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
  
  # output
  cat("Distance to", gene_type, "gene, filtering overlaping genes =", as.character(filter_overlap), "\n",
      "median =", median(abs(distance_to_nearest_elements$distance)), "\n", 
      "mean =", mean(abs(distance_to_nearest_elements$distance)), "\n")
  
  if(plot){
    
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
      ggsave(filename = paste0("MT2full_top100_distanceTo_", gene_type, "_genes_filterOverlap_", as.character(filter_overlap), ".pdf"), width = 20, height = 10)
    
  }
  
}

# calculates distance to nearest downstream genes (all or expressed), outputs mean and median + plot
distanceToDownstreamGene <- function(gene_type, plot){
  
  if(gene_type == "all"){
    genes_ranges <- knownGenes_gtf_gr
  }else{
    if(gene_type == "expressed"){
      genes_ranges <- knownGenes_expressed
    }
  }
  
  ######### distance to downstream genes
  ### elements which overlap genes
  # get elements 
  elements_overlaping_genes <- MT2OriginalGR[unique(queryHits(findOverlaps(MT2OriginalGR, genes_ranges, ignore.strand = T)))]
  
  # get distances (= 0)
  distance_to_overlaping_genes <- rep(0, length(elements_overlaping_genes))

  ### elements which are followed by gene downstream (gene downstream can be in any direction)
  # get elements 
  elements_followed_by_genes <- MT2OriginalGR[!(mcols(MT2OriginalGR)$repName %in% mcols(elements_overlaping_genes)$repName)]
  
  # get genes which follow elements
  genes_which_follow_elements <- genes_ranges[follow(x = elements_followed_by_genes, subject = genes_ranges, ignore.strand = T)]
  
  # get distances between elements and genes which follow them 
  distance_to_following_genes <- distance(elements_followed_by_genes, genes_which_follow_elements, ignore.strand = T)
  
  ###
  # put both distances to one vector
  distance_to_all <- c(distance_to_overlaping_genes, distance_to_following_genes)
  
  # calculate and print statistics
  cat("Distance to", gene_type, "genes downstream of MT2 full elements \n",
      "median =", (median(distance_to_all) / 1000), "kb\n", 
      "mean =", (mean(distance_to_all) / 1000), "kb\n")
  
  if(plot){
    
    # create data.frame
    distance_to_all_df <- data.frame(distance = distance_to_all, 
                                     gene_stream = "downstream")
    
    # plot histogram
    plot_hist <- 
      ggplot(data = distance_to_all_df, aes(x = distance)) +
      geom_histogram(binwidth = 1000)
    
    # boxplot downstream
    plot_box <- 
      ggplot(data = distance_to_all_df, aes(x = gene_stream, y = distance)) +
      geom_boxplot(fill = NA) +
      geom_jitter(color = "red", size = 1, height = 0.05) +
      coord_flip()
    
    # plot in grid
    plot_grid(plot_hist, plot_box, nrow = 2, align = "v") +
      ggsave(filename = paste0("MT2full_top100_distanceTo_", gene_type, "_downstream_genes.pdf"), width = 20, height = 10)
    
  }
  
  
}

distanceToDownstreamGene(gene_type = "all", plot = F)
distanceToDownstreamGene(gene_type = "expressed", plot = F)
