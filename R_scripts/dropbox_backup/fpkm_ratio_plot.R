setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/")

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

fpkm_all_filtered <- fpkm_all[!(fpkm_all[, "GV_KO_mean"] == 0 | fpkm_all[, "GV_WT_mean"] == 0 | fpkm_all[, "GV_WE"] == 0), ]
fpkm_all_filtered_no_zero <- fpkm_all_filtered[rowSums(fpkm_all_filtered == 0) == 0, ] 
fpkm_all_filtered_some_zero <- fpkm_all_filtered[rowSums(fpkm_all_filtered == 0) > 0, ] 

# colnames CNOT6: "GV_KO_mean"   "GV_WT_mean"   "MII_KO_mean"  "MII_WT_mean"  "OneC_KO_mean" "OneC_WT_mean" 
# colnames Fugaku: "OneCell_PA" "MII_PA" "OneCell_WE" "GV_WE" "MII_WE"

CNOT6_MII_KO_vs_MII_WT_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "MII_KO_mean", ratio2 = "MII_WT_mean")
CNOT6_MII_KO_vs_MII_WT_ratio_some_zero <- df_ratio(fpkm_all_filtered_some_zero, ratio1 = "MII_KO_mean", ratio2 = "MII_WT_mean")

CNOT6_1C_KO_vs_1C_WT_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "OneC_KO_mean", ratio2 = "OneC_WT_mean")
CNOT6_1C_KO_vs_1C_WT_ratio_some_zero <- df_ratio(fpkm_all_filtered_some_zero, ratio1 = "OneC_KO_mean", ratio2 = "OneC_WT_mean")

Fugaku_1C_PA_vs_MII_PA_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "OneCell_PA", ratio2 = "MII_PA")
Fugaku_1C_PA_vs_MII_PA_ratio_some_zero <- df_ratio(fpkm_all_filtered_some_zero, ratio1 = "OneCell_PA", ratio2 = "MII_PA")

Fugaku_1C_WE_vs_MII_WE_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "OneCell_WE", ratio2 = "MII_WE")
Fugaku_1C_WE_vs_MII_WE_ratio_some_zero <- df_ratio(fpkm_all_filtered_some_zero, ratio1 = "OneCell_WE", ratio2 = "MII_WE")

Fugaku_MII_WE_vs_GV_WE_ratio_no_zero <- df_ratio(fpkm_all_filtered_no_zero, ratio1 = "MII_WE", ratio2 = "GV_WE")
Fugaku_MII_WE_vs_GV_WE_ratio_some_zero <- df_ratio(fpkm_all_filtered_some_zero, ratio1 = "MII_WE", ratio2 = "GV_WE")


ratio_plot <- function(df_x, df_y, xlab, ylab, title){
  
  plot_df <- cbind(df_x, df_y)
  plot_df <- log2(plot_df)
  plot_df <- do.call(data.frame, lapply(plot_df, function(x) replace(x, is.infinite(x), NA))) #replaces Inf with NA
  colnames(plot_df) <- c("x", "y")
  
  current_plot <- ggplot() +
    geom_point(data = plot_df, size = 0.75, aes(x = x, y = y)) +
    scale_x_continuous(limits = c(-7.5, 7.5)) + 
    scale_y_continuous(limits = c(-7.5, 7.5)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title)
  
  return(current_plot)
}


plot1_no_zero <- ratio_plot(df_x = Fugaku_1C_PA_vs_MII_PA_ratio_no_zero, 
                            df_y = Fugaku_1C_WE_vs_MII_WE_ratio_no_zero, 
                            xlab = "1C/MII PA (log2(FPKM))", 
                            ylab = "1C/MII WE (log2(FPKM))", 
                            title = "1C/MII PA vs. 1C/MII WE")
plot1_some_zero <- ratio_plot(df_x = Fugaku_1C_PA_vs_MII_PA_ratio_some_zero, 
                            df_y = Fugaku_1C_WE_vs_MII_WE_ratio_some_zero, 
                            xlab = "1C/MII PA (log2(FPKM))", 
                            ylab = "1C/MII WE (log2(FPKM))", 
                            title = "1C/MII PA vs. 1C/MII WE")

plot2_no_zero <- ratio_plot(df_x = CNOT6_1C_KO_vs_1C_WT_ratio_no_zero, 
                            df_y = Fugaku_1C_WE_vs_MII_WE_ratio_no_zero, 
                            xlab = "1C KO/WT (log2(FPKM))", 
                            ylab = "1C/MII WE (log2(FPKM))", 
                            title = "1C KO/WT vs. 1C/MII WE")
plot2_some_zero <- ratio_plot(df_x = CNOT6_1C_KO_vs_1C_WT_ratio_some_zero, 
                            df_y = Fugaku_1C_WE_vs_MII_WE_ratio_some_zero, 
                            xlab = "1C KO/WT (log2(FPKM))", 
                            ylab = "1C/MII WE (log2(FPKM))", 
                            title = "1C KO/WT vs. 1C/MII WE")

plot3_no_zero <- ratio_plot(df_x = CNOT6_MII_KO_vs_MII_WT_ratio_no_zero, 
                            df_y = Fugaku_MII_WE_vs_GV_WE_ratio_no_zero, 
                            xlab = "MII KO/WT (log2(FPKM))", 
                            ylab = "MII/GV WE (log2(FPKM))", 
                            title = "MII KO/WT vs. MII/GV WE")
plot3_some_zero <- ratio_plot(df_x = CNOT6_MII_KO_vs_MII_WT_ratio_some_zero, 
                            df_y = Fugaku_MII_WE_vs_GV_WE_ratio_some_zero, 
                            xlab = "MII KO/WT (log2(FPKM))", 
                            ylab = "MII/GV WE (log2(FPKM))", 
                            title = "MII KO/WT vs. MII/GV WE")

plot <- list(plot1_no_zero, plot1_some_zero, plot2_no_zero, plot2_some_zero, plot3_no_zero, plot3_some_zero)

pdf("ratio_plots_fpkm_filtered.pdf")
invisible(lapply(plot, print))
dev.off()


order_table <- function(df){
  tx <- table(df)
  tx <- tx[order(tx, decreasing = T)]
  return(tx)
}

CNOT6_MII_KO_vs_MII_WT_ratio_no_zero_table <- head(order_table(CNOT6_MII_KO_vs_MII_WT_ratio_no_zero), n = 10)
CNOT6_1C_KO_vs_1C_WT_ratio_no_zero_table <- head(order_table(CNOT6_1C_KO_vs_1C_WT_ratio_no_zero), n = 10)
Fugaku_1C_PA_vs_MII_PA_ratio_no_zero_table <- head(order_table(Fugaku_1C_PA_vs_MII_PA_ratio_no_zero), n = 10)
Fugaku_1C_WE_vs_MII_WE_ratio_no_zero_table <- head(order_table(Fugaku_1C_WE_vs_MII_WE_ratio_no_zero), n = 10)
Fugaku_MII_WE_vs_GV_WE_ratio_no_zero_table <- head(order_table(Fugaku_MII_WE_vs_GV_WE_ratio_no_zero), n = 10)

# colnames CNOT6: "GV_KO_mean"   "GV_WT_mean"   "MII_KO_mean"  "MII_WT_mean"  "OneC_KO_mean" "OneC_WT_mean" 
# colnames Fugaku: "OneCell_PA" "MII_PA" "OneCell_WE" "GV_WE" "MII_WE"
a <- rownames(Fugaku_1C_PA_vs_MII_PA_ratio_no_zero)[which(as.character(Fugaku_1C_PA_vs_MII_PA_ratio_no_zero[, 1]) %in% names(Fugaku_1C_PA_vs_MII_PA_ratio_no_zero_table)[2])]
fpkm_all_filtered_no_zero[a, c("OneCell_PA", "MII_PA")]

counts_Fugaku <- read.csv("counts/counts_Fugaku.csv", row.names = 1)
counts_Fugaku[a, c("s_1cell.PA.bam", "s_MII.PA.bam")]
