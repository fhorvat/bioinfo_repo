library("ggplot2")
library("geneplotter")

byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}

dds <- DESeqDataSet(se[, 1:6], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[1:6]
dds_DESeq <- DESeq(dds)
res_GV <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_GV <- data.frame(res_GV)

dds <- DESeqDataSet(se[, 7:12], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[7:12]
dds_DESeq <- DESeq(dds)
res_MII <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_MII <- data.frame(res_MII)

dds <- DESeqDataSet(se[, 13:18], design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sampleTable_exp$ID[13:18]
dds_DESeq <- DESeq(dds)
res_1C <- results(dds_DESeq, contrast = c("Treatment.Control", "KO", "WT"))
res_1C <- data.frame(res_1C)

rpkm_old <- read.csv("rpkm_old.csv", row.names = 1)
colnames(rpkm_old) <- c("1C_PA", "MII_PA")

rpkm_new <- read.csv("rpkm_new.csv", row.names = 1)

rpkm_new_means <- byapply(rpkm_new, 3, rowMeans)
colnames(rpkm_new_means) <- c("GV_KO", "GV_WT", 
                              "MII_KO", "MII_WT", 
                              "1C_KO", "1C_WT")

rpkm_all <- merge(rpkm_new_means, rpkm_old, by = 0, all = T)
rownames(rpkm_all) <- rpkm_all$Row.names
rpkm_all <- rpkm_all[, -1]

rpkm_for_plot <- rpkm_all[, c(grep("1C", colnames(rpkm_all)), grep("MII", colnames(rpkm_all)))]
rpkm_for_plot[rpkm_for_plot == 0] <- NA
rpkm_for_plot$"x1C_MII_PA" <- rpkm_for_plot$"1C_PA" / rpkm_for_plot$"MII_PA"
rpkm_for_plot$"x1C_MII_WT" <- rpkm_for_plot$"1C_WT" / rpkm_for_plot$"MII_WT"
rpkm_for_plot$"x1C_MII_KO" <- rpkm_for_plot$"1C_KO" / rpkm_for_plot$"MII_KO"
rpkm_for_plot <- rpkm_for_plot[, 7:9]

rpkm_for_plot <- log2(rpkm_for_plot)
range_x <- range(c(na.omit(rpkm_for_plot$"x1C_MII_KO"), na.omit(rpkm_for_plot$"x1C_MII_PA")))
range_y <- range(c(na.omit(rpkm_for_plot$"x1C_MII_WT"), na.omit(rpkm_for_plot$"x1C_MII_PA")))

above <- 1  
genes_log2FC_above <- unique(unlist(lapply(list(res_1C, res_MII, res_GV), function(df) rownames(df)[which(abs(df$log2FoldChange) > above)])))
genes_log2FC_above <- rpkm_for_plot[which(rownames(rpkm_for_plot) %in% genes_log2FC_above), ]
genes_log2FC_above$gene_names <- rownames(genes_log2FC_above)

# RPKM 1C/MII PA vs. 1C/MII WT 
ggplot(rpkm_for_plot, aes(x = x1C_MII_PA, y = x1C_MII_WT)) +
  geom_point(size = 1) +
  geom_point(data = genes_log2FC_above, 
             aes(x = x1C_MII_PA, y = x1C_MII_WT), 
             color = "red3", 
             size = 3) +
  scale_x_continuous("log2 RKPM 1C/MII PA", limit = range_x) + 
  scale_y_continuous("log2 RPKM 1C/MII WT", limit = range_y) +
  ggtitle(paste("log2 RPKM 1C/MII PA vs. 1C/MII WT, genes with |log2FC| >", above, "in KO (GV, MII, 1C) in red"))
savepng(paste0("2D_log2_RPKM_1C_MII_PA_vs_1C_MII_WT_log2FC_above_", above), width = 1000)

# RPKM 1C/MII KO vs. 1C/MII WT 
ggplot(rpkm_for_plot, aes(x = x1C_MII_KO, y = x1C_MII_WT)) +
  geom_point(size = 1) +
  geom_point(data = genes_log2FC_above, 
             aes(x = x1C_MII_KO, y = x1C_MII_WT), 
             color = "red3", 
             size = 3) +
  scale_x_continuous("log2 RKPM 1C/MII KO", limit = range_x) + 
  scale_y_continuous("log2 RPKM 1C/MII WT", limit = range_y) +
  ggtitle(paste("log2 RPKM 1C/MII KO vs. 1C/MII WT, genes with |log2FC| >", above, "in KO (GV, MII, 1C) in red"))
savepng(paste0("2D_log2_RPKM_1C_MII_KO_vs_1C_MII_WT_log2FC_above_", above), width = 1000)

# RPKM 1C/MII KO vs. 1C/MII PA 
ggplot(rpkm_for_plot, aes(x = x1C_MII_KO, y = x1C_MII_PA)) +
  geom_point(size = 1) +
  geom_point(data = genes_log2FC_above, 
             aes(x = x1C_MII_KO, y = x1C_MII_PA), 
             color = "red3", 
             size = 3) +
  scale_x_continuous("log2 RKPM 1C/MII KO", limit = range_x) + 
  scale_y_continuous("log2 RPKM 1C/MII PA", limit = range_y) +
  ggtitle(paste("log2 RPKM 1C/MII KO vs. 1C/MII PA, genes with |log2FC| >", above, "in KO (GV, MII, 1C) in red"))
savepng(paste0("2D_log2_RPKM_1C_MII_KO_vs_1C_MII_PA_log2FC_above_", above), width = 1000)
