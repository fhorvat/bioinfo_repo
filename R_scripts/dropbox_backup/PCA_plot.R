library("ggplot2")
setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/DESeq_analysis_Svoboda/Prague/results_polyA_seq/counts_rpkm/2cell_Deng_2014")

all_exons_gene_counts <- read.csv("all_exons_counts.csv", stringsAsFactors = F, row.names = 1)  

plot_PCA <- function(object){
  
  # filteres and orders data by highest variance in first 500 rows 
  ntop <- 500
  rv <- genefilter::rowVars(object) 
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  # calculates pca
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  # makes data.frame
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample_name = colnames(object))
  attr(d, "percentVar") <- percentVar[1:2]
  
  # plots
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "sample_name")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
}

plot_PCA(all_exons_gene_counts)

