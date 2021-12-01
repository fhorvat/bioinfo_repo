library("ggplot2")
library("geneplotter")

rpkm_new <- read.csv("rpkm_new.csv", row.names = 1)
sampleTable <- read.csv("/common/WORK/fhorvat/RNA_Seq_CNOT6L/CNOT6L_sample list_11919R_2015-10-29.csv", header = T)
sampleTable <- sampleTable[1:18, ]

rpkm_new[rpkm_new == 0] <- NA
rpkm_new <- do.call(data.frame, lapply(rpkm_new, function(x) replace(x, is.infinite(x), NA)))
rownames(rpkm_new) <- rownames(rpkm_new)
rpkm_new <- rpkm_new[rowSums(is.na(rpkm_new)) == 0, ]
rpkm_new <- log2(rpkm_new)
colnames(rpkm_new) <- c(colnames(rpkm_new)[1:12], paste0("oneC_KO", 1:3), paste0("oneC_WT", 1:3))

pca <- prcomp(t(rpkm_new), scale = F)
scores <- data.frame(colnames(rpkm_new), 
                     pca$x[,1:2], 
                     TimeCourse = sampleTable$Time.Course, 
                     TreatmentControl = sampleTable$Treatment.Control, 
                     TimeCourse_Treatment = paste(sampleTable$Time.Course, sampleTable$Treatment.Control, sep = "_"))

ggplot(scores, aes(PC1, PC2, color = TreatmentControl, shape = TimeCourse)) + 
  geom_point(size = 3)
savepng("PCA_log2_RPKM_all_samples2", width = 1000, asp = 0.75)
