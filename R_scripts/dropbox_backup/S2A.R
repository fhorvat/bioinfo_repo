rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## FUNCTIONS
insertFreqByGenomeElement <- function(category){
  
  # filter by category
  rpmsk_element <- rptmsk_viz_filtered_insert[rptmsk_viz_filtered_insert$category_unique == category]
  
  # find overlaps
  overlaps_promoter <- findOverlaps(rpmsk_element, knownGenes_promoters)
  overlaps_exon <- findOverlaps(rpmsk_element, knownGenes_exons)
  overlaps_intron <- findOverlaps(rpmsk_element, knownGenes_introns)
  overlaps_intergenic <- findOverlaps(rpmsk_element, knownGenes_intergenic)
  
  # filter overlaps: promoter > exon > intron > intergenic
  overlaps_intergenic <- overlaps_intergenic[-which(queryHits(overlaps_intergenic) %in% c(queryHits(overlaps_promoter),
                                                                                          queryHits(overlaps_exon),
                                                                                          queryHits(overlaps_intron)))]
  overlaps_intron <- overlaps_intron[-which(queryHits(overlaps_intron) %in% c(queryHits(overlaps_promoter),
                                                                              queryHits(overlaps_exon)))]
  overlaps_exon <- overlaps_exon[-which(queryHits(overlaps_exon) %in% queryHits(overlaps_promoter))]
  
  # count overlaps by regions
  region_count <- 
    lapply(X = list(overlaps_intergenic, 
                    overlaps_promoter, 
                    overlaps_intron, 
                    overlaps_exon), 
           FUN = function(X) length(unique(queryHits(X)))) %>% 
    set_names(., c("intergenic", "promoters", "introns", "exons")) %>% 
    unlist()
  
  # insert frequency by total genomic element length in kb 
  insert_freq <- region_count / region_length
  
  return(insert_freq)
  
}

######################################################## READ DATA
### knownGenes table
knownGenes_gtf <- makeTxDbFromGFF("UCSC_knownGene_mm10_20161126.gtf.gz")

# genes
knownGenes_genes <- genes(knownGenes_gtf, single.strand.genes.only = FALSE)
strand(knownGenes_genes) <- "*"

# intergenic regions
knownGenes_intergenic <- gaps(unlist(knownGenes_genes))

# promoter regions (1kb upstream)
knownGenes_promoters <- unlist(knownGenes_genes)
end(knownGenes_promoters) <- start(knownGenes_promoters) - 1
start(knownGenes_promoters) <- start(knownGenes_promoters) - 1000

# introns
knownGenes_introns <- intronsByTranscript(knownGenes_gtf)
knownGenes_introns <- unlist(knownGenes_introns)
strand(knownGenes_introns) <- "*"

# exons
knownGenes_exons <- exonsBy(knownGenes_gtf)
knownGenes_exons <- unlist(knownGenes_exons)
strand(knownGenes_exons) <- "*"

# length of all genomic regions in kb
region_length <- 
  lapply(X = list(knownGenes_intergenic, knownGenes_promoters, knownGenes_introns, knownGenes_exons), FUN = function(X) sum(width(reduce(X)) / 1000)) %>% 
  set_names(., c("intergenic", "promoters", "introns", "exons")) %>% 
  unlist()

### repeatMaskerVIZ table
rptmsk_viz <- 
  read_delim("UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz", delim = "\t") %>% 
  dplyr::select(seqnames = genoName, start = genoStart, end = genoEnd, strand, repName, viz_id = id) 

# table with all LTRs, create patterns for matching in repeatMasker table
LTR_repName_pattern <- 
  read_csv("LTR_families_data.csv") %$% 
  repName %>% 
  paste0("^", .) %>% 
  paste(., collapse = "|")

# table with category for each element (categories = MLT, ORR1, MT, MT2 or other)
LTR_families <- read_csv("LTR_families_data_expanded.csv") 

######################################################## MAIN CODE
# coordinates of all LTRs from repeatMaskerVIZ table
rptmsk_viz_filtered <- 
  rptmsk_viz[str_detect(rptmsk_viz$repName, LTR_repName_pattern), ] %>% 
  dplyr::filter(!str_detect(repName, "^LTR16A2$|^LTR41C$|^LTR75_1$"))

# split -int from LTRs
rptmsk_viz_filtered_int <- rptmsk_viz_filtered[str_detect("-int", rptmsk_viz_filtered$repName), ]
rptmsk_viz_filtered_ltr <- rptmsk_viz_filtered[!str_detect("-int", rptmsk_viz_filtered$repName), ]

# combine inserts with at least one LTR with solo LTRs
rptmsk_viz_ltr_int_id <- 
  left_join(rptmsk_viz_filtered_int, rptmsk_viz_filtered_ltr, by = "viz_id") %>% 
  dplyr::select(seqnames.x, start.x, end.x, repName.x, start.y, end.y, repName.y, viz_id) %>% 
  dplyr::filter(complete.cases(.)) %$%
  viz_id %>% 
  c(., rptmsk_viz_filtered_ltr$viz_id) %>% 
  unique()

rptmsk_viz_filtered %<>%  
  dplyr::filter(viz_id %in% rptmsk_viz_ltr_int_id)

# inserts with all sub-elements
rptmsk_viz_filtered_ID <- 
  rptmsk_viz_filtered %>%
  left_join(LTR_families, by = "repName") %>% 
  group_by(viz_id) %>%
  dplyr::select(repName, category, viz_id) %>%
  summarise(repName_unique = str_c(repName, collapse = "|"), 
            category_paste = str_c(category, collapse = "|")) %>% 
  mutate(category_unique = 
           category_paste %>% 
           str_split(., "\\|") %>% 
           sapply(X = ., FUN = function(X) str_c(unique(X), collapse = "|"))) %>% 
  dplyr::select(viz_id, repName_unique, category_unique)

# coordinates of whole inserts
rptmsk_viz_filtered_insert <- 
  makeGRangesFromDataFrame(rptmsk_viz_filtered, keep.extra.columns = T) %>% 
  split(mcols(.)$viz_id) %>% 
  range(ignore.strand = T) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  mutate(viz_id = as.integer(rownames(.))) %>% 
  left_join(rptmsk_viz_filtered_ID, by = "viz_id") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# insert frequency in distinctive genomic regions for each category of LTR (MLT, ORR1, MT, MT2, other)
insert_frequency <- 
  lapply(unique(rptmsk_viz_filtered_insert$category_unique), insertFreqByGenomeElement) %>% 
  do.call(rbind, .) %>% 
  set_rownames(., unique(rptmsk_viz_filtered_insert$category_unique)) %>% 
  as.data.frame() %>% 
  mutate(category = factor(rownames(.), levels = c("MLT", "ORR1", "MT", "MT2", "other"))) %>% 
  arrange(category) %>% 
  dplyr::select(c(5, 1:4)) %T>% 
  write_csv(x = ., path = "LTR_insert_freq_in_genomic_regions.csv")