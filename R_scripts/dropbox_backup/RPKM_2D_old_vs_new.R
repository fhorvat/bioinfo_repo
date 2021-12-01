library("ggplot2")
library("geneplotter")

rpkm_old <- read.csv("rpkm_old.csv", row.names = 1)
colnames(rpkm_old) <- c("oneC_PA", "MII_PA")
rpkm_new <- read.csv("rpkm_new.csv", row.names = 1)
rpkm_new_means <- byapply(rpkm_new, 3, rowMeans)
colnames(rpkm_new_means) <- c("GV_KO", "GV_WT", 
                              "MII_KO", "MII_WT", 
                              "oneC_KO", "oneC_WT")

rpkm_all <- merge(rpkm_new_means, rpkm_old, by = 0, all = T)
rownames(rpkm_all) <- rpkm_all$Row.names
rpkm_all <- rpkm_all[, -1]

rpkm_for_plot <- rpkm_all
rpkm_for_plot[rpkm_for_plot == 0] <- NA
rpkm_for_plot <- rpkm_for_plot[rowSums(is.na(rpkm_for_plot)) < 1, ]
rpkm_for_plot <- log2(rpkm_for_plot)
range_x <- range(rpkm_for_plot$oneC_WT, rpkm_for_plot$MII_WT)
range_y <- range(rpkm_for_plot$oneC_PA, rpkm_for_plot$MII_PA)

fit <- lm(oneC_WT ~ 0 + oneC_PA, data = rpkm_for_plot)
lb <- paste("R^2 == ", format(summary(fit)$adj.r.squared, digits = 4))   
ggplot(rpkm_for_plot, aes(x = oneC_WT, y = oneC_PA)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = coef(fit)[1], color = "red") +
  annotate("text", x = range_x[2] - 1, y = range_y[1], label = lb, parse = TRUE) +
  scale_x_continuous("log2 RKPM 1C WT", limit = range_x) + 
  scale_y_continuous("log2 RPKM 1C PA", limit = range_y) +
  ggtitle("")
savepng("log2_RKPM_1C_WT_vs_1C_PA", width = 1000)


fit <- lm(MII_WT ~ 0 + MII_PA, data = rpkm_for_plot)
lb <- paste("R^2 == ", format(summary(fit)$adj.r.squared, digits = 4))   
ggplot(rpkm_for_plot, aes(x = MII_WT, y = MII_PA)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = coef(fit)[1], color = "red") +
  annotate("text", x = range_x[2] - 1, y = range_y[1], label = lb, parse = TRUE) +
  scale_x_continuous("log2 RKPM MII WT", limit = range_x) + 
  scale_y_continuous("log2 RPKM MII PA", limit = range_y) +
  ggtitle("")
savepng("log2_RKPM_MII_WT_vs_MII_PA", width = 1000)
