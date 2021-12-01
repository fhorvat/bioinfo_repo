setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/")

# CNOT6L
log2fc_MII_KO_vs_MII_WT_ratio <- read.csv("diff_exp/DESeq2_KOvsWT_MII.csv", row.names = 1)
log2fc_MII_KO_vs_MII_WT_ratio <- data.frame(log2fc_MII_KO_vs_MII_WT_ratio[, "log2FoldChange"], 
                                            row.names = rownames(log2fc_MII_KO_vs_MII_WT_ratio))
colnames(log2fc_MII_KO_vs_MII_WT_ratio) <- "log2fc_MII_KO_vs_MII_WT"

log2fc_1C_KO_vs_1C_WT_ratio <- read.csv("diff_exp/DESeq2_KOvsWT_1C.csv", row.names = 1)
log2fc_1C_KO_vs_1C_WT_ratio <- data.frame(log2fc_1C_KO_vs_1C_WT_ratio[, "log2FoldChange"], 
                                          row.names = rownames(log2fc_1C_KO_vs_1C_WT_ratio))
colnames(log2fc_1C_KO_vs_1C_WT_ratio) <- "log2fc_1C_KO_vs_1C_WT"


# FUGAKU
log2fc_1C_PA_vs_MII_PA_ratio <- read.csv("diff_exp/DESeq2_1CvsMII_PA.csv", row.names = 1)
log2fc_1C_PA_vs_MII_PA_ratio <- data.frame(log2fc_1C_PA_vs_MII_PA_ratio[, "log2FoldChange"], 
                                           row.names = rownames(log2fc_1C_PA_vs_MII_PA_ratio))
colnames(log2fc_1C_PA_vs_MII_PA_ratio) <- "log2fc_1C_PA_vs_MII_PA"

log2fc_1C_WE_vs_MII_WE_ratio <- read.csv("diff_exp/DESeq2_1CvsMII_WE.csv", row.names = 1)
log2fc_1C_WE_vs_MII_WE_ratio <- data.frame(log2fc_1C_WE_vs_MII_WE_ratio[, "log2FoldChange"], 
                                           row.names = rownames(log2fc_1C_WE_vs_MII_WE_ratio))
colnames(log2fc_1C_WE_vs_MII_WE_ratio) <- "log2fc_1C_WE_vs_MII_WE"

log2fc_MII_WE_vs_GV_WE_ratio <- read.csv("diff_exp/DESeq2_MIIvsGV_WE.csv", row.names = 1)
log2fc_MII_WE_vs_GV_WE_ratio <- data.frame(log2fc_MII_WE_vs_GV_WE_ratio[, "log2FoldChange"], 
                                           row.names = rownames(log2fc_MII_WE_vs_GV_WE_ratio))
colnames(log2fc_MII_WE_vs_GV_WE_ratio) <- "log2fc_MII_WE_vs_GV_WE"


ratio_plot <- function(df_x, df_y, xlab, ylab, title, xlim, ylim){
  
  plot_df <- cbind(df_x, df_y)
  # plot_df <- do.call(data.frame, lapply(plot_df, function(x) replace(x, is.infinite(x), NA))) #replaces Inf with NA
  colnames(plot_df) <- c("x", "y")
  
  current_plot <- ggplot() +
    geom_point(data = plot_df, size = 0.75, aes(x = x, y = y)) +
    scale_x_continuous(limits = xlim) + 
    scale_y_continuous(limits = ylim) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title)
  
  return(current_plot)
}

ratio_plot(log2fc_1C_KO_vs_1C_WT_ratio, 
           log2fc_1C_WE_vs_MII_WE_ratio, 
           xlab = "1C KO/WT (log2FC)", 
           ylab = "1C/MII WE (log2FC)", 
           title = "1C KO/WT vs. 1C/MII WE", 
           xlim = c(-10, 10),
           ylim = c(-5, 20))

ratio_plot(log2fc_1C_PA_vs_MII_PA_ratio, 
           log2fc_1C_WE_vs_MII_WE_ratio,
           xlab = "1C/MII PA (log2FC)", 
           ylab = "1C/MII WE (log2FC)", 
           title = "1C/MII PA vs. 1C/MII WE", 
           xlim = c(-10, 20),
           ylim = c(-5, 20))

ratio_plot(log2fc_MII_KO_vs_MII_WT_ratio, 
           log2fc_MII_WE_vs_GV_WE_ratio,
           xlab = "MII KO/WT (log2FC)", 
           ylab = "MII/GV WE (log2FC)", 
           title = "MII KO/WT vs. MII/GV WE", 
           xlim = c(-10, 10),
           ylim = c(-5, 20))



