library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("geneplotter")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("genefilter")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("gage")
library("pathview")
library("dplyr")
library("magrittr")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Analysis/diffExp")

# exons by genes
ebg <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
exons_width <- sapply(width(ebg), sum)

# files lnc Nov2016 experiment
sample_table_lnc2016_full <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/RNAseq_2016_11_23_sampleTable.csv", header = T) %>%
  dplyr::select(ID = sample, Time.Course = Sample, Treatment.Control = Type) %>%
  mutate(Treatment.Control = gsub(" B6", "", Treatment.Control), 
         ID = gsub("_16.*", "", ID), 
         name = paste(gsub(".*_[1-8]_", "", ID), Treatment.Control, sep = "_"), 
         experiment = "lnc_Nov2016") %>% 
  left_join(data.frame(sample_paths = list.files("/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                                 pattern = "*.bam$",
                                                 recursive = T,
                                                 full.names = T), 
                       logs_paths = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/lnc_2016_RNAseq/Data/Mapped/STAR_mm10", 
                                               pattern = "*Log.final.out", 
                                               recursive = T, 
                                               full.names = T), 
                       stringsAsFactors = F) %>%
              mutate(ID = gsub("^/.*/|_16.*", "", sample_paths)) %>%
              rowwise() %>%
              mutate(lib_size = (as.integer(read.delim(logs_paths, header = F, stringsAsFactors = F)[8, 2]) / 10^6)), 
            by = "ID")

################################################################## counting reads overlaping with exons from knownGene 
bamfiles_lnc2016 <- BamFileList(sample_table_lnc2016_full[, "sample_paths"], yieldSize = 2000000)
register(MulticoreParam())
se_lnc <- summarizeOverlaps(features = ebg, 
                            reads = bamfiles_lnc2016, 
                            mode = "Union", 
                            singleEnd = TRUE, 
                            ignore.strand = TRUE)

# filter
lnc_sample <- "Lnc5"
sample_table_filter <- sample_table_lnc2016_full[grepl(lnc_sample, sample_table_lnc2016_full$Treatment.Control) | sample_table_lnc2016_full$Treatment.Control == "WT", ]
# sample_table_filter$Treatment.Control2 <- ifelse(grepl("Lnc", sample_table_lnc2016_full$Treatment.Control), "Lnc", "WT")
sample_table_filter <- sample_table_filter[c(2, 3), ]
sample_index_filter <- as.numeric(rownames(sample_table_filter))
se_filter <- se_lnc[, sample_index_filter]
colData(se_filter) <- DataFrame(sample_table_filter)

################################################################## fpkm table
fpkm_df <- as.data.frame(assay(se_filter))
colnames(fpkm_df) <- sample_table_filter$"ID"
fpkm_df$width <- exons_width
invisible(lapply(X = as.character(sample_table_filter$"ID"), 
                 FUN = function(X) fpkm_df[, X] <<- fpkm_df[, X] / (sample_table_filter[sample_table_filter$ID == X, "lib_size"] * (fpkm_df$width / 1000))))
fpkm_df <- fpkm_df[, -ncol(fpkm_df)]
colnames(fpkm_df) <- paste(sample_table_filter[, "name"], rep(1:2, (nrow(sample_table_filter) / 2)), "rpkm", sep = "_")
fpkm_df$entrezID <- rownames(fpkm_df)

################################################################## 
dds <- DESeqDataSet(se_filter, design = ~Treatment.Control)
dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- colnames(fpkm_df)
rld <- rlog(dds, blind = F)

################################################################## results
dds_DESeq <- DESeq(dds)
res <- results(dds_DESeq, contrast = c("Treatment.Control", lnc_sample, "WT"))
res <- res[complete.cases(res), ]

# join with fpkm
res <- as.data.frame(res)
res$entrezID <- rownames(res)
res <- left_join(res, fpkm_df, by = "entrezID")

# assigning gene names and gene symbols
res$gene_symbol <- mapIds(org.Mm.eg.db, keys = res$entrezID, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
res$gene_name <- mapIds(org.Mm.eg.db, keys = res$entrezID, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")
res <- res[, c("entrezID", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "BC6_Lnc5_1_rpkm",  "BC7_WT_2_rpkm", "gene_name")]

res_up <- res[order(res$log2FoldChange, decreasing = T), ]
res_down <- res[order(res$log2FoldChange, decreasing = F), ]

# getting and writing significantly differentaly expressed genes
# res_signf <- subset(res, padj < 0.1)
# write.csv(res_signf, paste0("diffExp_", lnc_sample, "_vs_WT_signf.csv")) 
write.csv(res_up, "diffExp_upregulated_Lnc5BC6_vs_WTBC7.csv") 
write.csv(res_down, "diffExp_downregulated_Lnc5BC6_vs_WTBC7.csv") 

# MA plot
res_plot <- data.frame(mean = rowMeans(res[, c("BC6_Lnc5_1_rpkm",  "BC7_WT_2_rpkm")]), 
                       lfc = res$log2FoldChange,
                       sig = "BLCK",
                       gene = res$gene_symbol,
                       stringsAsFactors = F)
res_plot[res_plot$mean > 1000 & res_plot$mean < 1500, "sig"] <- "RD"
res_plot[res_plot$sig != "RD", "gene"] <- NA
res_plot$mean_log <- ifelse(res_plot$mean > 1000, log10(res_plot$mean), res_plot$mean)

MA_plot <- ggplot(data = res_plot, aes(x = mean_log, y = lfc, color = sig, alpha = sig)) + 
  geom_point(size = 0.1) +
  geom_abline(slope = 0, color = "red3", alpha = 0.7) +
  geom_text(aes(label = gene), vjust = "inward", hjust = "inward", nudge_y = -0.2, nudge_x = 0.2, size = 3, color = "black") +
  scale_y_continuous(limits = c(-4, 4)) + 
  scale_colour_manual(values = c(BLCK = "gray32", RD = "red3", BLU = "blue3")) +
  scale_alpha_manual(values = c(BLCK = 0.3, RD = 1, BLU = 1)) +
  guides(color = FALSE, alpha = FALSE) +
  xlab("mean expression") + 
  ylab("log2FoldChange lnc5/WT") +
  ggsave("MAplot_Lnc5BC6_vs_WTBC7_8.png")

################################################################## plots
# mannualy set colors to the stages
color_pallete <- hue_pal()(8)[c(1, 5, 7, 6, 3, 2)]
names(color_pallete) <- c("WT", "KO", "Lnc5", "Lnc12", "Lnc21", "Lnc")

plotPCAandDistHeatmap <- function(plot_data){
  
  ## data for plot (fpkm or rld)
  if(plot_data == "fpkm"){
    data_df <- fpkm_df
  }else{
    if(plot_data == "rld"){
      data_df <- assay(rld)
    }else{
      stop("wrong data")
    }
  }
  
  #### PCA
  # filteres and orders data by highest variance in first 500 rows 
  ntop <- 500
  rv <- genefilter::rowVars(data_df) 
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  # calculates pca
  pca <- prcomp(t(data_df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  # makes data.frame
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  name = sample_table_filter$"name", 
                  Treatment.Control = sample_table_filter$"Treatment.Control", 
                  lib_size = sample_table_filter$"lib_size", 
                  Time.Course = sample_table_filter$"Time.Course")
  d$Treatment.Control <- factor(d$Treatment.Control2, levels = c("WT", "KO", "Lnc5", "Lnc12", "Lnc21", "Lnc"))
  
  # plot
  # show_col(hue_pal()(8))
  ggplot(d, aes(PC1, PC2)) + 
    geom_point(aes(size = lib_size, color = Treatment.Control)) +
    geom_text(aes(label = name), vjust = "inward", hjust = "inward", nudge_y = -0.2, nudge_x = 0.2, size = 2, color = "black") +
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    scale_color_manual(values = color_pallete) + 
    ggsave(paste0("PCA_", plot_data, "_lnc2016_", lnc_sample, "_WT.png"))
  
  ### Distance matrix (heatmap)
  sampleDists <- dist(t(data_df))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  png(paste0("dist_", plot_data, "_lnc2016_", lnc_sample, "_WT.png"), width = 1000, height = 1000)
  pheatmap(sampleDistMatrix, 
           clustering_distance_rows = sampleDists, 
           clustering_distance_cols = sampleDists, 
           col = colors)
  dev.off()
  
  return("done ploting")
}

### MA plot
png(paste0("MAplot_lnc2016_", lnc_sample, "_WT.png"), width = 1000, height = 1000)
plotMA(res, ylim = c(-2, 2), cex = 1)
dev.off()

### PCA and distance heatmap
invisible(lapply(list("fpkm", "rld"), plotPCAandDistHeatmap))

################################################################### pathview
# # KEGG set for mouse
# kg.mouse <- kegg.gsets("mouse")
# kegg.gs <- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
# 
# # samples
# sample_list <- "lnc5_vs_WT"
# 
# # gage
# deseq2.fc <- res$log2FoldChange
# names(deseq2.fc) <- rownames(res)
# fc.kegg.p <- gage(deseq2.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
# 
# # upregulated pathways
# sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
# path.ids <- rownames(fc.kegg.p$greater)[sel]
# path.ids <- substr(path.ids, 1, 8)
# df_upregulated <- data.frame(q.value = fc.kegg.p$greater[sel, "q.val"])
# if (nrow(df_upregulated) > 1){
#   df_upregulated$number <- 1 : nrow(df_upregulated)  
# }
# write.csv(df_upregulated, paste0("./gage_pathview/", "gage_", sample_list, "_upreg.csv"))
# pv_upregulated <- sapply(path.ids, function(pid) pathview(gene.data = deseq2.fc, 
#                                                           pathway.id = pid, 
#                                                           species = "mmu",
#                                                           limit = list(gene = 3, cpd = 2),
#                                                           bins = list(gene = 20, cpd = 20),
#                                                           out.suffix = paste0(sample_list, "_upreg"), 
#                                                           kegg.dir = "./gage_pathview"))
# 
# # downregulated pathways
# sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
# path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
# path.ids.l <- substr(path.ids.l, 1, 8)
# df_downregulated <- data.frame(q.value = fc.kegg.p$less[sel.l, "q.val"])
# if (nrow(df_downregulated) > 1){
#   df_downregulated$number <- 1 : nrow(df_downregulated)  
# }
# write.csv(df_downregulated, paste0("./gage_pathview/", "gage_", sample_list, "_upreg.csv"))
# pv_downregulated <- sapply(path.ids.l, function(pid) pathview(gene.data = deseq2.fc, 
#                                                               pathway.id = pid, 
#                                                               species = "mmu",
#                                                               limit = list(gene = 3, cpd = 2),
#                                                               bins = list(gene = 20, cpd = 20),
#                                                               out.suffix = paste0(sample_list, "_downreg"), 
#                                                               kegg.dir = "./gage_pathview/"))
# 
