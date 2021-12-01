library(GenomicRanges)
library(rtracklayer)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/L1_elements/mapping/references")
reads <- import.bed(con = "L1_all_20160418.bed.gz")

# writes LINE-1 elements longer than 6000 bp into BED file
reads_filtered <- reads[width(reads) > 6000]
reads_df <- as.data.frame(reads_filtered)
reads_df <- reads_df[, -4]
reads_df <- reads_df[, c("seqnames", "start", "end", "name", "score", "strand")]
write.table(reads_df, "L1_6000bp.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# writes only first 500 bp-s of LINE-1 elements longer than 6000 bp into BED file
reads_filtered_first_500bp <- reads_filtered
plus_strand <- IRanges(start = start(reads_filtered_first_500bp), width = 500)
minus_strand <- IRanges(end = end(reads_filtered_first_500bp), width = 500)

for (X in 1:length(reads_filtered_first_500bp)){
  if(as.logical(strand(reads_filtered_first_500bp[X]) == "+")){
    ranges(reads_filtered_first_500bp[X]) = plus_strand[X]
  } else{
    ranges(reads_filtered_first_500bp[X]) = minus_strand[X]
  }  
}
reads_filtered_first_500bp_df <- as.data.frame(reads_filtered_first_500bp)
reads_filtered_first_500bp_df <- reads_filtered_first_500bp_df[, -4]
reads_filtered_first_500bp_df <- reads_filtered_first_500bp_df[, c("seqnames", "start", "end", "name", "score", "strand")]
write.table(reads_filtered_first_500bp_df, "L1_6000bp_500bp_5UTR.bed", quote = F, sep = "\t", row.names = F, col.names = F)

# 12 column .BED
reads_filtered_12col <- reads_filtered
mcols(reads_filtered_12col)$thickStart <- start(reads_filtered_12col)
mcols(reads_filtered_12col)$thickEnd <- end(reads_filtered_12col)
mcols(reads_filtered_12col)$itemRgb <- 1
mcols(reads_filtered_12col)$blockCount <- 1
mcols(reads_filtered_12col)$blockSizes <- width(reads_filtered_12col)
mcols(reads_filtered_12col)$blockStarts <- 0

reads_filtered_12col_df <- as.data.frame(reads_filtered_12col)
reads_filtered_12col_df <- reads_filtered_12col_df[, c("seqnames", "start", "end", "name", 
                                                       "score", "strand", "thickStart", "thickEnd",
                                                       "itemRgb", "blockCount", "blockSizes", "blockStarts")]
write.table(reads_filtered_12col_df, "L1_6000bp_12column.bed", quote = F, sep = "\t", row.names = F, col.names = F)


reads_filtered_first_500bp_12col <- reads_filtered_first_500bp
mcols(reads_filtered_first_500bp_12col)$thickStart <- start(reads_filtered_first_500bp_12col)
mcols(reads_filtered_first_500bp_12col)$thickEnd <- end(reads_filtered_first_500bp_12col)
mcols(reads_filtered_first_500bp_12col)$itemRgb <- 1
mcols(reads_filtered_first_500bp_12col)$blockCount <- 1
mcols(reads_filtered_first_500bp_12col)$blockSizes <- width(reads_filtered_first_500bp_12col)
mcols(reads_filtered_first_500bp_12col)$blockStarts <- 0

reads_filtered_first_500bp_12col_df <- as.data.frame(reads_filtered_first_500bp_12col)
reads_filtered_first_500bp_12col_df <- reads_filtered_first_500bp_12col_df[, c("seqnames", "start", "end", "name",
                                                                               "score", "strand", "thickStart", "thickEnd",
                                                                               "itemRgb", "blockCount", "blockSizes", "blockStarts")]
write.table(reads_filtered_first_500bp_12col_df, "L1_6000bp_500bp_5UTR_12column.bed", quote = F, sep = "\t", row.names = F, col.names = F)
