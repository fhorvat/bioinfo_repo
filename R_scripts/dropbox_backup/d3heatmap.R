library("RColorBrewer")
library("d3heatmap")

rm(list = ls())

# server
samples_fpkm <- read.table("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/oocyte_specific_genes/Fugaku_knownGene_fpkm.txt", sep = "\t", stringsAsFactors = F, row.names = NULL)[, -1]

# local
setwd("C:/Users/fhorvat/Dropbox/Praksa bioinfo/Projekti/test")
samples_fpkm <- read.table("Fugaku_knownGene_fpkm.txt", sep = "\t", stringsAsFactors = F, row.names = NULL)[, -1]

# sample names
samples <- c("s_GV.WE", "s_MII.WE", "s_1cell.WE", "s_2cell.WE", "s_4cell.WE", "s_Molura.WE", "s_Blast.WE")

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

samples_fpkm_filtered <- samples_fpkm_filtered[, c(1, 10, 2:3, 11, 4:6, 12, 7:8, 13)]
colnames(samples_fpkm_filtered)[3:12] <- paste0(colnames(samples_fpkm_filtered)[3:12], "_fpkm")
rownames(samples_fpkm_filtered) <- samples_fpkm_filtered$symbol

########################################################################## heatmap
# order by GV + MII + 1C
samples_fpkm_filtered_heatmap <- samples_fpkm_filtered[, c("gene_id", paste0(samples, "_fpkm"))]
samples_fpkm_filtered_heatmap <- samples_fpkm_filtered_heatmap[order(rowSums(samples_fpkm_filtered[, c("s_GV.WE_fpkm", "s_MII.WE_fpkm", "s_1cell.WE_fpkm")]), decreasing = T), ]

# log(FPKM + 2)
samples_fpkm_filtered_heatmap <- samples_fpkm_filtered_heatmap[, paste0(samples, "_fpkm")]
samples_fpkm_filtered_heatmap <- log2(samples_fpkm_filtered_heatmap + 2)
samples_fpkm_filtered_heatmap <- t(samples_fpkm_filtered_heatmap)
rownames(samples_fpkm_filtered_heatmap) <- c("GV", "MII", "1-cell", "2-cell", "4-cell",  "morula", "blastocyst")

# plot only 100 samples
samples_fpkm_filtered_heatmap_sample <- samples_fpkm_filtered_heatmap[, 1:100]

# color pallete
bk <- seq(range(samples_fpkm_filtered_heatmap_sample)[1], range(samples_fpkm_filtered_heatmap_sample)[2], length = 50)
hmcols <- rev(colorRampPalette(c("red", "yellow", "black"))(length(bk) - 1))

# plot
d3heatmap(samples_fpkm_filtered_heatmap_sample,
          Rowv = F,
          Colv = F, 
          main = "log(FPKM + 2)",
          colors = hmcols,
          show_grid = F, 
          xaxis_font_size = 8)
