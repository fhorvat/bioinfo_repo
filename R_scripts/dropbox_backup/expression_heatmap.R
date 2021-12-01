library("DESeq2")
library("GenomicRanges")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
library("VennDiagram")
library("geneplotter")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

rm(list = ls())
options(bitmapType = 'cairo')

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/oocyte_specific_genes")

######################################################################### getting count of reads over knownGene table
# samples paths
samples <- c("s_GV.WE", "s_MII.WE", "s_1cell.WE", "s_2cell.WE", "s_4cell.WE", "s_Molura.WE", "s_Blast.WE")
samples_path <- "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/"
filenames <- paste0(samples_path, samples, "/", samples, ".bam")
filenames_logs <- paste0(samples_path, samples, "/", samples, "Log.final.out")

# # making TxDb object from knownGene from UCSC
# ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
# gene_name <- mapIds(org.Mm.eg.db,
#                     keys = names(ebg),
#                     column ="GENENAME",
#                     keytype = "ENTREZID",
#                     multiVals = "first")
# 
# # filtering non-protein coding genes
# gene_name <- names(grep("non-protein", gene_name, value = T))
# ebg <- ebg[!(names(ebg) %in% gene_name)]
# 
# # counting overlaps
# bamfiles <- BamFileList(filenames[3], yieldSize = 2000000)
# se <- summarizeOverlaps(features = ebg, 
#                         reads = bamfiles, 
#                         mode = "Union", 
#                         singleEnd = FALSE, 
#                         ignore.strand = TRUE)
# 
# # getting data.frame of counts
# samples_counts <- as.data.frame(assay(se))
# colnames(samples_counts) <- samples
# samples_counts$gene_id <- names(ebg)
# samples_counts$width <- sapply(width(ebg), sum)
# rownames(samples_counts) <- NULL
# write.table(samples_counts, "Fugaku_knownGene_counts.txt", sep = "\t", col.names = T, quote = F)

samples_counts <- read.delim("Fugaku_knownGene_counts.txt", stringsAsFactors = F)
                    
########################################################################## calculating FPKM 
# getting library size in millions of reads
number_of_reads <- sapply(X = 1:length(filenames_logs), function(X) as.integer(read.delim(filenames_logs[X], header = F, stringsAsFactors = F)[8, 2]))
names(number_of_reads) <- samples
number_of_reads <- number_of_reads / 10^6

# calculating FPKM from counts
samples_fpkm <- samples_counts
for (x in 1:length(number_of_reads)){
  samples_fpkm[, samples[x]] <- samples_fpkm[, samples[x]] / (number_of_reads[samples[x]] * (samples_fpkm$width / 1000))
}

# samples_fpkm$symbol <- mapIds(org.Mm.eg.db, 
#                               keys = rownames(samples_fpkm), 
#                               column ="SYMBOL",
#                               keytype = "ENTREZID",
#                               multiVals = "first")
# write.table(samples_fpkm, "Fugaku_knownGene_fpkm.txt", sep = "\t", col.names = T, quote = F)

# calculating average FPKM (maternal/zygotic/embrional)
samples_fpkm$maternal <-  rowMeans(samples_fpkm[, samples[1:2]]) 
samples_fpkm$zygotic <- rowMeans(samples_fpkm[, samples[4:5]])
samples_fpkm$embrional <- rowMeans(samples_fpkm[, samples[6:7]])

########################################################################## filtering
# keep genes which are:
# - max. in maternal
# - min. in embrional
# - embrional < 5% maternal
# - at least 1 FPKM in GV or MII
# - less than 1 FPKM in blastocyst

# getting max and min values in maternal, zygotic and embrional FPKM
samples_fpkm$max <- apply(samples_fpkm[, c("maternal", "zygotic", "embrional")], 1, which.max)
samples_fpkm$min <- apply(samples_fpkm[, c("maternal", "zygotic", "embrional")], 1, which.min)

# filtering based on criteria above
samples_fpkm_filtered <- samples_fpkm[samples_fpkm$max == 1 &
                                        samples_fpkm$min == 3 &
                                        samples_fpkm$embrional < (0.05 * samples_fpkm$maternal) & 
                                        (samples_fpkm$s_GV.WE > 1 | samples_fpkm$s_MII.WE > 1) & 
                                        samples_fpkm$s_Blast.WE < 1, ] 

########################################################################## gene list, different gene names, writting table
samples_fpkm_filtered$gene_symbol <- mapIds(org.Mm.eg.db,
                                       keys = as.character(samples_fpkm_filtered$gene_id),
                                       column = "SYMBOL",
                                       keytype = "ENTREZID",
                                       multiVals = "first")
samples_fpkm_filtered$gene_name <- mapIds(org.Mm.eg.db,
                                          keys = as.character(samples_fpkm_filtered$gene_id),
                                          column = "GENENAME",
                                          keytype = "ENTREZID",
                                          multiVals = "first")

samples_fpkm_filtered <- samples_fpkm_filtered[, c(1, 15, 2:3, 10, 4:6, 11, 7:8, 12, 16)]
colnames(samples_fpkm_filtered)[3:12] <- paste0(colnames(samples_fpkm_filtered)[3:12], "_fpkm")
# write.csv(samples_fpkm_filtered, "gene_list_with_fpkm_2.csv", row.names = F)

########################################################################## heatmap
# order by GV + MII + 1C
samples_fpkm_filtered_heatmap <- samples_fpkm_filtered[, c("gene_id", paste0(samples, "_fpkm"))]
samples_fpkm_filtered_heatmap <- samples_fpkm_filtered_heatmap[order(rowSums(samples_fpkm_filtered[, c("s_GV.WE_fpkm", "s_MII.WE_fpkm", "s_1cell.WE_fpkm")]), decreasing = T), ]

# log(FPKM + 2)
samples_fpkm_filtered_heatmap <- samples_fpkm_filtered_heatmap[, paste0(samples, "_fpkm")]
samples_fpkm_filtered_heatmap <- log2(samples_fpkm_filtered_heatmap + 2)
samples_fpkm_filtered_heatmap <- t(samples_fpkm_filtered_heatmap)
rownames(samples_fpkm_filtered_heatmap) <- c("GV", "MII", "1-cell", "2-cell", "4-cell",  "morula", "blastocyst")
  
# color pallete
bk <- seq(range(samples_fpkm_filtered_heatmap)[1], range(samples_fpkm_filtered_heatmap)[2], length = 50)
hmcols <- rev(colorRampPalette(c("red", "yellow", "black"))(length(bk) - 1))

# plot and save
png("pheatmap_fpkm_2.png", width = 1400, height = 700)
print(pheatmap(samples_fpkm_filtered_heatmap,
               main = "log(FPKM + 2)",
               col = hmcols, 
               breaks = bk,
               cluster_row = F, 
               cluster_cols = F,
               show_rownames = T, 
               show_colnames = F))
dev.off()



############################################################1############## plotting Venn diagram of MII and 1C with FPKM > 1 
# Venn diagram plot three samples
s_MII_gene_id <- samples_fpkm[samples_fpkm$"s_MII.WE" > 1, "gene_id"]
s_2C_gene_id <- samples_fpkm[samples_fpkm$"s_2cell.WE" > 1, "gene_id"]

grid.newpage()
venn.plot <- venn.diagram(list(samples_fpkm$gene_id, s_MII_gene_id, s_2C_gene_id), 
                          NULL, 
                          fill = NULL, 
                          cex = 2, 
                          cat.fontface = 4, 
                          cat.cex = 1.7,
                          category.names = c("protein coding genes", "MII", "2C"))
grid.draw(venn.plot)
savepng("MII_2C_fpkm_above_1", width = 1000, asp = 1)
