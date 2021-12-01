library("ggplot2")
library("scales")
library("org.Mm.eg.db")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MousePolyA/2cell/Hamazaki_2015")
tx_fpkm_original <- read.csv("FPMK_tx_width.csv", stringsAsFactors = F, row.names = 1)

tx_fpkm <- tx_fpkm_original
tx_fpkm <- tx_fpkm[rowSums(tx_fpkm[ ,names(tx_fpkm)[grep("FPKM", names(tx_fpkm))]]) > 0, ]
tx_fpkm <- tx_fpkm[, c(1, 4:ncol(tx_fpkm))]
tx_fpkm$average_fpkm <- rowMeans(tx_fpkm[, 2:ncol(tx_fpkm)])


# plot with highlight
tx_fpkm <- tx_fpkm_original
tx_fpkm <- tx_fpkm[rowSums(tx_fpkm[ ,names(tx_fpkm)[grep("FPKM", names(tx_fpkm))]]) > 0, ]
tx_fpkm <- tx_fpkm[, c(1, 4:ncol(tx_fpkm))]
tx_fpkm$average_fpkm <- rowMeans(tx_fpkm[, 2:ncol(tx_fpkm)])
tx_fpkm$highlight <- ifelse(tx_fpkm$average_fpkm > 2000, "highlight", "normal")
textdf <- tx_fpkm[tx_fpkm$highlight == "highlight", ]
mycolours <- c("highlight" = "red", "normal" = "grey50")

ggplot(tx_fpkm, aes(x = width_sum, y = average_fpkm, fill = highlight)) +
  geom_jitter(width = 1.5, size = 3, pch = 21, colour = "White") +
  geom_text(data = textdf, aes(x = width_sum + 20000, y = average_fpkm, label = rownames(textdf))) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c()) + 
  scale_y_continuous(limits = c())
  

# plot with all samples
tx_fpkm <- tx_fpkm_original
tx_fpkm <- tx_fpkm[, c(1, 4:ncol(tx_fpkm))]
# tx_fpkm <- tx_fpkm[sample(x = nrow(tx_fpkm), size = 100, replace = F), ]
df_plot <- data.frame(fpkm = unlist(lapply(X = 2:ncol(tx_fpkm), function(X) tx_fpkm[, X])),
                      width = tx_fpkm$width_sum,
                      sample_name = unlist(lapply(colnames(tx_fpkm)[2:ncol(tx_fpkm)], rep, nrow(tx_fpkm))))

ggplot(df_plot, aes(width, fpkm, fill = sample_name)) + 
  geom_jitter(width = 0, size = 3, pch = 21, colour = "White") +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(limits = c()) + 
  scale_y_continuous(limits = c())


# density plot of samples
# on lobsang
tx_fpkm <- tx_fpkm_original
tx_fpkm <- tx_fpkm[rowSums(tx_fpkm[ ,names(tx_fpkm)[grep("FPKM", names(tx_fpkm))]]) > 0, ]
tx_fpkm <- tx_fpkm[, c(1, 4:ncol(tx_fpkm))]

df_plot <- data.frame(fpkm = unlist(lapply(X = 2:ncol(tx_fpkm), function(X) tx_fpkm[, X])),
                      sample_name = unlist(lapply(colnames(tx_fpkm)[2:ncol(tx_fpkm)], rep, nrow(tx_fpkm))))

write.csv(df_plot, "fpkm_for_density_plot.csv", row.names = T)

# local
setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/DESeq_analysis_Svoboda/Prague/results_polyA_seq")
df_plot <- read.csv("fpkm_for_density_plot.csv", row.names = 1)
df_plot <- df_plot[df_plot$fpkm > 0, ]

ggplot(df_plot) + 
  geom_density(aes(x = log10(fpkm), 
                   group = sample_name, 
                   color = sample_name, 
                   fill = sample_name), 
               alpha = I(1/3))
