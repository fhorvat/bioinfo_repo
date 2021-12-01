## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite("chromPlot")

library("chromPlot")
library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
library("tidyverse")
library("data.table")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/annotation")

################################################################## getting genes, exons, introns, promoters and intergenic regions from knownGene tables
# knownGenes table
knownGenes_gtf <- makeTxDbFromGFF("/common/WORK/fhorvat/reference/mm10/UCSC_knownGenes_mm10_20160513.gtf.gz")

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
rptmsk_viz <- read.delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMaskerVIZ_OutBaseline_20160824.txt.gz", stringsAsFactors = F, header = T)
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
  select(element_name, VIZ_ID) %>%
  summarise(element_name_paste = paste(element_name, collapse = "|"))

# add full element names to ranges
rptmsk_viz_filtered_insert_df <- as.data.frame(rptmsk_viz_filtered_insert)
rptmsk_viz_filtered_insert_df$VIZ_ID <- as.integer(rptmsk_viz_filtered_insert_df$VIZ_ID)
rptmsk_viz_filtered_insert_df <- left_join(rptmsk_viz_filtered_insert_df, rptmsk_viz_filtered_ID, by = "VIZ_ID")
rptmsk_viz_filtered_insert <- makeGRangesFromDataFrame(rptmsk_viz_filtered_insert_df, keep.extra.columns = T)

