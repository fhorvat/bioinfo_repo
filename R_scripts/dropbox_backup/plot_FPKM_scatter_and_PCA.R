library("ggplot2")
library("scales")
library("org.Mm.eg.db")

plot_topN_genes <- function(gene_counts, N, plot_boxplot){
  gene_counts <- gene_counts[order(rowMeans(gene_counts), decreasing = T), ] 
  gene_counts_top_100 <- gene_counts[1:N, ]
  
  gene_names <- mapIds(org.Mm.eg.db,
                       keys = rownames(gene_counts_top_100),
                       column ="SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

  df_plot <- data.frame(unlist(lapply(X = 1:ncol(gene_counts), function(X) gene_counts_top_100[, X])),
                        gene_names, 
                        unlist(lapply(colnames(gene_counts), rep, N)))
  colnames(df_plot) <- c("count", "gene", "sample_name")
  
  if(plot_boxplot == T){
    ggplot(df_plot, aes(gene, count, fill = sample_name)) + 
      geom_jitter(width = 1.5, size = 3, pch = 21, colour = "White") +
      geom_boxplot() +
      scale_y_continuous(labels = comma) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
      guides(fill = guide_legend(override.aes = list(size = 5)))  +
      scale_x_discrete(limits = gene_names) 
  } else {
    ggplot(df_plot, aes(gene, count, fill = sample_name)) + 
      geom_jitter(width = 1.5, size = 3, pch = 21, colour = "White") +
      scale_y_continuous(labels = comma) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
      guides(fill = guide_legend(override.aes = list(size = 5)))  +
      scale_x_discrete(limits = gene_names) 
  }
}

plot_PCA <- function(object){
  ntop <- 500
  rv <- genefilter::rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample_name = colnames(object))
  attr(d, "percentVar") <- percentVar[1:2]
  
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "sample_name")) + 
    geom_point(size = 3) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
}

plot_and_save <- function(df_for_plot){
  current_df <- read.csv(file = df_for_plot, stringsAsFactors = F, row.names = 1)
  
  # scatter plot
  scatter_width <- 19.1
  scatter_height <- 7.61
  plot_name <- paste0("plot_", 
                      gsub("_counts.csv|_rpkm.csv", "", df_for_plot), 
                      "_scatter_", 
                      gsub("^.*_exons_|.csv", "", df_for_plot),
                      "_top_100_genes.pdf")
  current_plot <- plot_topN_genes(current_df, N = 100, plot_boxplot = F)
  ggsave(plot_name, plot = current_plot, width = scatter_width, height =  scatter_height)
  
  #PCA plot
  PCA_width <- 9.03
  PCA_height <- 7.61
  plot_name <- paste0("plot_", 
                      gsub("_counts.csv|_rpkm.csv", "", df_for_plot), 
                      "_PCA_", 
                      gsub("^.*_exons_|.csv", "", df_for_plot),
                      ".pdf")
  current_plot <- plot_PCA(current_df)
  ggsave(plot_name, plot = current_plot, width = PCA_width, height =  PCA_height)
  return("plotting done")
}

plot_in_dir <- function(working_directory){
  setwd(working_directory)
  sapply(c("all_exons_counts.csv", 
           "all_exons_rpkm.csv",  
           "last_exons_counts.csv", 
           "last_exons_rpkm.csv"),
         plot_and_save)  
}

work_dir <- c("/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/2cell/Hamazaki_2015/",
              "/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/MII/Hamazaki_2015/")
sapply(work_dir, plot_in_dir)
