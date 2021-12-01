library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
library("dplyr")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/annotation")

################################################################## getting genes, exons, introns, promoters and intergenic regions from knownGene tables
# knownGenes table
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz")

#get genes
knownGenes_genes <- genes(knownGenes_gtf, single.strand.genes.only = FALSE)
strand(knownGenes_genes) <- "*"

# get intergenic regions
knownGenes_intergenic <- gaps(unlist(knownGenes_genes))

# get promoter regions (1kb upstream)
knownGenes_promoters <- unlist(knownGenes_genes)
end(knownGenes_promoters) <- start(knownGenes_promoters) - 1
start(knownGenes_promoters) <- start(knownGenes_promoters) - 1000

#get introns
knownGenes_introns <- intronsByTranscript(knownGenes_gtf)
knownGenes_introns <- unlist(knownGenes_introns)
strand(knownGenes_introns) <- "*"

#get exons
knownGenes_exons <- exonsBy(knownGenes_gtf)
knownGenes_exons <- unlist(knownGenes_exons)
strand(knownGenes_exons) <- "*"


################################################################## get MT2/ORR/MLT LTRs and inserts from repeatMasker
rptmsk_viz <- read.delim("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz", stringsAsFactors = F, header = T)
rptmsk_viz <- rptmsk_viz[, c("genoName", "genoStart", "genoEnd", "strand", "repName", "id")]
colnames(rptmsk_viz) <- c("seqnames", "start", "end", "strand", "element_name", "VIZ_ID")

# filter for MT, ORR and MLT
rptmsk_viz_filtered <- rptmsk_viz[grep("MT|ORR|MLT", rptmsk_viz$element_name), ]
rptmsk_viz_filtered <- rptmsk_viz_filtered[-grep("MMTV-int|MLTR", rptmsk_viz_filtered$element_name), ]

################################################################## split int from LTRs, check if whole insert with -int also contains at least one LTR 
# split int from LTRs
rptmsk_viz_filtered_int <- rptmsk_viz_filtered[grep("-int", rptmsk_viz_filtered$element_name), ]
rptmsk_viz_filtered_ltr <- rptmsk_viz_filtered[-grep("-int", rptmsk_viz_filtered$element_name), ]

# get inserts with at least one LTR
rptmsk_viz_filtered_int_with_ltr <- left_join(rptmsk_viz_filtered_int, rptmsk_viz_filtered_ltr, by = "VIZ_ID")
rptmsk_viz_filtered_int_with_ltr <- rptmsk_viz_filtered_int_with_ltr[, c("seqnames.x", "start.x", "end.x", "element_name.x", "start.y", "end.y", "element_name.y", "VIZ_ID")]
rptmsk_viz_filtered_int_with_ltr <- rptmsk_viz_filtered_int_with_ltr[complete.cases(rptmsk_viz_filtered_int_with_ltr), ]

# take VIZ_ID of those inserts with at least one LTR and combine them with VIZ_ID of LTRs
rptmsk_viz_ltr_int_id <- unique(c(rptmsk_viz_filtered_ltr$VIZ_ID, rptmsk_viz_filtered_int_with_ltr$VIZ_ID))

################################################################## get ranges of whole insert based on VIZ id, filter, get full names of inserts
# get ranges of whole insert based on VIZ id
rptmsk_viz_filtered_insert <- makeGRangesFromDataFrame(rptmsk_viz_filtered, keep.extra.columns = T)
rptmsk_viz_filtered_insert <- split(rptmsk_viz_filtered_insert, mcols(rptmsk_viz_filtered_insert)$VIZ_ID)
rptmsk_viz_filtered_insert <- range(rptmsk_viz_filtered_insert, ignore.strand = T)
rptmsk_viz_filtered_insert <- unlist(rptmsk_viz_filtered_insert)
rptmsk_viz_filtered_insert$VIZ_ID <- names(rptmsk_viz_filtered_insert)

# filter ranges of inserts without LTRs (only -int)
rptmsk_viz_filtered_insert <- rptmsk_viz_filtered_insert[names(rptmsk_viz_filtered_insert) %in% rptmsk_viz_ltr_int_id]
names(rptmsk_viz_filtered_insert) <- NULL

# get names of elements (LTRs and ints) which make one insert
rptmsk_viz_filtered_ID <- 
  rptmsk_viz_filtered %>%
  group_by(VIZ_ID) %>%
  dplyr::select(element_name, VIZ_ID) %>%
  summarise(element_name_paste = paste(element_name, collapse = "|"))

# add full element names to ranges
rptmsk_viz_filtered_insert_df <- as.data.frame(rptmsk_viz_filtered_insert)
rptmsk_viz_filtered_insert_df$VIZ_ID <- as.integer(rptmsk_viz_filtered_insert_df$VIZ_ID)
rptmsk_viz_filtered_insert_df <- left_join(rptmsk_viz_filtered_insert_df, rptmsk_viz_filtered_ID, by = "VIZ_ID")
rptmsk_viz_filtered_insert <- makeGRangesFromDataFrame(rptmsk_viz_filtered_insert_df, keep.extra.columns = T)


################################################################## get frequencies of insertions into genome elements

insertFreqByGenomeElement <- function(repeat_id){
  
  rpmsk_element <- rptmsk_viz_filtered_insert[grep(repeat_id, rptmsk_viz_filtered_insert$element_name_paste)]
  if (repeat_id == "MT"){
    rpmsk_element <- rpmsk_element[-grep("MT2", rpmsk_element$element_name_paste)]
  }
  
  # finding overlaps
  overlaps_promoter <- findOverlaps(rpmsk_element, knownGenes_promoters)
  overlaps_exon <- findOverlaps(rpmsk_element, knownGenes_exons)
  overlaps_intron <- findOverlaps(rpmsk_element, knownGenes_introns)
  overlaps_intergenic <- findOverlaps(rpmsk_element, knownGenes_intergenic)
  
  # filtering overlaps in order promoter > exon > intron > intergenic
  overlaps_intergenic <- overlaps_intergenic[-which(queryHits(overlaps_intergenic) %in% c(queryHits(overlaps_promoter),
                                                                                          queryHits(overlaps_exon),
                                                                                          queryHits(overlaps_intron)))]
  overlaps_intron <- overlaps_intron[-which(queryHits(overlaps_intron) %in% c(queryHits(overlaps_promoter),
                                                                              queryHits(overlaps_exon)))]
  overlaps_exon <- overlaps_exon[-which(queryHits(overlaps_exon) %in% queryHits(overlaps_promoter))]
  
  
  # counting overlaps by regions
  region_count <- lapply(X = list(overlaps_intergenic, overlaps_promoter, overlaps_intron, overlaps_exon), 
                         FUN = function(X) length(unique(queryHits(X))))
  names(region_count) <- c("intergenic", "promoters", "introns", "exons")
  region_count <- unlist(region_count)
  
  # getting insert frequency by total genomic element length in kb 
  insert_freq <- region_count / region_length
  return(insert_freq)
}

# get length of all genomic regions in kb
region_length <- lapply(X = list(knownGenes_intergenic, knownGenes_promoters, knownGenes_introns, knownGenes_exons), FUN = function(X) sum(width(reduce(X)) / 1000))
names(region_length) <- c("intergenic", "promoters", "introns", "exons")
region_length <- unlist(region_length)

insert_frequency <- lapply(c("ORR", "MT", "MT2", "MLT"), insertFreqByGenomeElement)
insert_frequency <- do.call(rbind, insert_frequency)
rownames(insert_frequency) <- c("ORR", "MT", "MT2", "MLT")
write.csv(x = insert_frequency, file = "MT_MT2_ORR_MLT_insertFreqByGenomeElements.csv")

