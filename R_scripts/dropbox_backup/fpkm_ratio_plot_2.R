setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/")
setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/Svoboda/CNOT6L/rpkm_ratio_with_Fugaku")

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

df_ratio <- function(df, ratio1, ratio2){
  df <- data.frame(df[, ratio1] / df[, ratio2], row.names = rownames(df))
  colnames(df) <- paste(ratio1, ratio2, "ratio", sep = "_")
  return(df)
}

fpkm_MII <- read.csv("FPKM/fpkm_MII.csv", row.names = 1)
fpkm_1C <- read.csv("FPKM/fpkm_1C.csv", row.names = 1)
fpkm_GV <- read.csv("FPKM/fpkm_GV.csv", row.names = 1)
fpkm_CNOT6 <- cbind(fpkm_GV, fpkm_MII, fpkm_1C)
fpkm_CNOT6 <- byapply(fpkm_CNOT6, 3, rowMeans)
colnames(fpkm_CNOT6) <- c("GV_KO_mean", "GV_WT_mean", "MII_KO_mean", "MII_WT_mean", "OneC_KO_mean", "OneC_WT_mean")

fpkm_Fugaku <- read.csv("FPKM/fpkm_Fugaku.csv", row.names = 1)
colnames(fpkm_Fugaku) <- c("OneCell_PA", "MII_PA", "OneCell_WE", "GV_WE", "MII_WE")

fpkm_all <- cbind(fpkm_CNOT6, fpkm_Fugaku)
fpkm_all["226777", ]
# fpkm of GV in CNOT6 WT&KO > 4
# fpkm_all_filtered <- fpkm_all[fpkm_all[, "GV_KO_mean"] > 4 & fpkm_all[, "GV_WT_mean"] > 4, ]

# fpkm of GV in Fugaku GV_WE > 4
fpkm_all_filtered <- fpkm_all[fpkm_all[, "GV_WE"] > 4, ]

fpkm_all_filtered <- fpkm_all

# remove those with fpkm = 0 in any stage
fpkm_all_filtered_no_zero <- fpkm_all_filtered[rowSums(fpkm_all_filtered == 0) == 0, ] 

# colnames CNOT6: "GV_KO_mean"   "GV_WT_mean"   "MII_KO_mean"  "MII_WT_mean"  "OneC_KO_mean" "OneC_WT_mean" 
# colnames Fugaku: "OneCell_PA" "MII_PA" "OneCell_WE" "GV_WE" "MII_WE"

CNOT6_1C_KO_vs_1C_WT_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "OneC_KO_mean", ratio2 = "OneC_WT_mean")
Fugaku_1C_PA_vs_MII_PA_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "OneCell_PA", ratio2 = "MII_PA")
Fugaku_1C_WE_vs_MII_WE_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "OneCell_WE", ratio2 = "MII_WE")


plot1_df <- cbind(Fugaku_1C_PA_vs_MII_PA_ratio_no_zero, Fugaku_1C_WE_vs_MII_WE_ratio_no_zero)
plot1_df <- log2(plot1_df)
plot1_df_names <- rownames(plot1_df[plot1_df[, 1] <= -1 & plot1_df[, 2] <= -1, ])

# red = all genes with:
# - 1C/MII PA ratio <= 1
# - 1C/MII WE ratio <= 1 
plot1_df_red <- plot1_df[plot1_df_names, ]

# blue = genes which are blue in plot 2 (= have ratio CNOT6L 1C KO/WT > 1)
plot1_df_blue <- plot1_df[plot2_df_blue_names, ]
colnames(plot1_df) <- c("x", "y")
colnames(plot1_df_red) <- c("x", "y")
colnames(plot1_df_blue) <- c("x", "y")

current_plot1 <- ggplot() +
  geom_point(data = plot1_df, size = 0.75, aes(x = x, y = y)) +
  geom_point(data = plot1_df_red, color = "red", size = 0.75, aes(x = x, y = y)) +
  geom_point(data = plot1_df_blue, color = "blue", size = 0.75, aes(x = x, y = y)) +
  scale_x_continuous(limits = c(-7.5, 7.5)) + 
  scale_y_continuous(limits = c(-7.5, 7.5)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab("1C/MII PA (log2(FPKM))") +
  ylab("1C/MII WE (log2(FPKM))") +
  ggtitle("1C/MII PA vs. 1C/MII WE")  

plot2_df <- cbind(CNOT6_1C_KO_vs_1C_WT_ratio_no_zero, Fugaku_1C_WE_vs_MII_WE_ratio_no_zero)
plot2_df <- log2(plot2_df)

# red = all genes which are red in plot 1
plot2_df_red <- plot2_df[plot1_df_names, ]

# blue = all red from plot 1 with ratio CNOT6L 1C KO/WT > 1
plot2_df_blue_names <- rownames(plot2_df_red[plot2_df_red[, 1] >= 1, ])
plot2_df_blue <- plot2_df_red[plot2_df_blue_names, ]
colnames(plot2_df) <- c("x", "y")
colnames(plot2_df_red) <- c("x", "y")
colnames(plot2_df_blue) <- c("x", "y")

current_plot2 <- ggplot() +
  geom_point(data = plot2_df, size = 0.75, aes(x = x, y = y)) +
  geom_point(data = plot2_df_red, color = "red", size = 0.75, aes(x = x, y = y)) +
  geom_point(data = plot2_df_blue, color = "blue", size = 0.75, aes(x = x, y = y)) +
  scale_x_continuous(limits = c(-7.5, 7.5)) + 
  scale_y_continuous(limits = c(-7.5, 7.5)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab("1C KO/WT (log2(FPKM))") +
  ylab("1C/MII WE (log2(FPKM))") +
  ggtitle("1C KO/WT vs. 1C/MII WE")  

plot <- list(current_plot1, current_plot2)

pdf("ratio_plots_fpkm_filtered_GV_fpkm_more_than_4_in_Fugaku_labeled_red_blue.pdf")
invisible(lapply(plot, print))
dev.off()

ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm9.knownGene, by = "gene")

gene_list <- rownames(plot1_df_red)
gene_list_ranges <- data.frame(range(ebg[gene_list]))
gene_list_ranges <- gene_list_ranges[, c("group_name", "seqnames", "start", "end")]
gene_list_ranges$symbol <- mapIds(org.Mm.eg.db, 
                                  keys = gene_list_ranges$group_name,
                                  column ="SYMBOL",
                                  keytype = "ENTREZID",
                                  multiVals = "first")
gene_list_ranges$gene_name <- mapIds(org.Mm.eg.db, 
                                     keys = gene_list_ranges$group_name,
                                     column ="GENENAME",
                                     keytype = "ENTREZID",
                                     multiVals = "first")

gene_list_pos <- plot1_df_red[gene_list_ranges$group_name, ]
gene_list_table <- cbind(gene_list_ranges, gene_list_pos)
colnames(gene_list_table) <- c("entrezID", "chr", "start", "end", "symbol", "gene_name", "x_pos", "y_pos")
gene_list_table <- gene_list_table[order(gene_list_table$y_pos, decreasing = T), ]
rownames(gene_list_table) <- NULL

write.csv(gene_list_table, "gene_list_red.csv", row.names = F)