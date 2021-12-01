#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 5. 4. 2017.  
### AUTHOR: Filip Horvat
### NOTE:
###### lnc5 comparison
# parts of code used to compare 2 lnc5 samples (so no replicates)
###### lnc5 comparison

rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp")

################################################################################### LIBRARIES
################################################################################### 
library(dplyr)
library(stringr)
library(readr)
library(magrittr)
library(ggplot2)
library(reshape2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(AnnotationDbi)
library(goseq)

library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(GO.db)

################################################################################### PATH VARIABLES
################################################################################### 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/Mapped/STAR_mm10"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp/results"

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/documentation/lncKO_2016_library_size.txt"

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

################################################################################### SOURCE FILES
################################################################################### 
source(file.path(lib_path, "GffToGRanges.R"))

################################################################################### FUNCTIONS
################################################################################### 

################################################################################### SCRIPT PARAMS
################################################################################### 

################################################################################### TABLES
################################################################################### 
# experiment table
sample_table <- 
  read_delim(sample_table_path, delim = "\t") %>% 
  dplyr::mutate(treatment = str_replace(name, "BC.*_", ""))

# unique treatment
treatments <- unique(sample_table$treatment)

# read in the ENSEMBL gene table, take protein coding genes
ensembl_gtf <- read_delim(file = mm10_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 

################################################################################### MAIN CODE
################################################################################### 
### take transcript with most exons for each gene
# get exons for each transcript
gtf_trans <- GffToGRanges(ensembl_gtf, "exon")

# get unique gene_id/transcript_id combinations
gids <- unique(values(gtf_trans)[c('gene_id', 'transcript_id')])

# splits exon ranges on transcripts
gtf_trans <- split(gtf_trans, gtf_trans$transcript_id)

# orders transcripts based on number of exons in transcripts
gtf_trans <- gtf_trans[order(elementNROWS(gtf_trans), decreasing = T)] 

# keeps only first transcript of each gene (the one with most exons)
gtf_trans <- gtf_trans[!duplicated(gids$gene_id[match(names(gtf_trans), gids$transcript_id)])]

# get total length of all exons for each transcript
trans_width <- 
  sum(width(gtf_trans)) %>% 
  tibble(transcript_id = names(.), width = .)
  
# get total length of all exons for each transcript
trans_width <- 
  gtf_trans %>%  
  unlist(.) %>% 
  as_data_frame() %>% 
  dplyr::select(gene_id, transcript_id) %>% 
  unique(.) %>% 
  left_join(sum(width(gtf_trans)) %>% 
              tibble(transcript_id = names(.), width = .))

################################################################################### 
### get count of reads 
# counting overlaps
bamfiles <- Rsamtools::BamFileList(sample_table$sample_path, yieldSize = 2000000)
register(MulticoreParam())
se <- GenomicAlignments::summarizeOverlaps(features = gtf_trans,
                                           reads = bamfiles,
                                           mode = "Union",
                                           singleEnd = TRUE,
                                           ignore.strand = TRUE)

# RPKM normalization
rpkm_df <- 
  as.data.frame(assay(se)) %>% 
  magrittr::set_colnames(sample_table$name) %>% 
  tibble::rownames_to_column(var = "transcript_id") %>% 
  tidyr::gather(key = name, value = counts, -transcript_id) %>% 
  dplyr::left_join(., sample_table[, c("name", "library_size")], by = "name") %>% 
  dplyr::left_join(., trans_width, by = "transcript_id") %>% 
  dplyr::mutate(library_size = round(library_size / 1E6, 2), 
                width = round(width / 1E3, 2), 
                rpm = (counts / library_size), 
                rpkm = (rpm / width)) %>% 
  dplyr::select(transcript_id, name, rpkm) %>% 
  tidyr::spread(key = name, value = rpkm)

# mean RPKM in treatment 
rpkm_mean <- 
  tidyr::gather(data = rpkm_df, key = name, value = rpkm, -transcript_id) %>% 
  dplyr::left_join(., sample_table[, c("name", "treatment")], by = "name") %>% 
  dplyr::group_by(transcript_id, treatment) %>% 
  dplyr::summarise(rpkm = round(mean(rpkm), 2)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = treatment, value = rpkm)
  
# set column data
colData(se) <- DataFrame(sample_table)

# DESeq2 differential expression - KO vs. WT
dds <- DESeqDataSet(se, design = ~treatment)

###### lnc5 comparison
dds <- DESeqDataSet(se, design = ~name)
###### lnc5 comparison

dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sample_table$ID

###### lnc5 comparison
dds <- dds[, c("HT3LFBGXY_5_BC5", "HT3LFBGXY_6_BC6")]
dds$name <- droplevels(dds$name)
###### lnc5 comparison

dds <- DESeq(dds)
# saveRDS(dds, file = "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp/results/dds_design_name.RDS")

###### lnc5 comparison
res_all <- 
  results(dds, contrast = c("name", "BC5_Lnc5", "BC6_Lnc5")) %>% 
  tibble::as_tibble(x = .) %>% 
  tibble::rownames_to_column(var = "transcript_id") %>% 
  dplyr::arrange(log2FoldChange) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::left_join(., rpkm_df[, c("transcript_id", "BC5_Lnc5", "BC6_Lnc5")], by = "transcript_id") %>% 
  data.table::setnames(., c("BC5_Lnc5", "BC6_Lnc5"), c("BC5_Lnc5_RPKM", "BC6_Lnc5_RPKM")) %>% 
  dplyr::mutate(gene_symbol = mapIds(EnsDb.Mmusculus.v79, keys = transcript_id, column = "SYMBOL", keytype = "TXID", multiVals = "first"),
                gene_biotype = mapIds(EnsDb.Mmusculus.v79, keys = transcript_id, column = "GENEBIOTYPE", keytype = "TXID", multiVals = "first"))

# add gene name
gene_name <- mapIds(org.Mm.eg.db, keys = res_all$gene_symbol, column = "GENENAME", keytype = "SYMBOL", multiVals = "first")
gene_name[sapply(gene_name, is.null)] <- NA
gene_name <- unname(unlist(gene_name))

# add GO ID
GO_ID <- mapIds(org.Mm.eg.db, keys = res_all$gene_symbol, column = "GO", keytype = "SYMBOL", multiVals = "first") 
GO_ID[sapply(GO_ID, is.null)] <- NA
GO_ID <- unname(unlist(GO_ID))

# add GO term
GO_term <- mapIds(GO.db, keys = GO_ID, column = "TERM", keytype = "GOID", multiVals = "first")
GO_term[sapply(GO_term, is.null)] <- NA
GO_term <- unname(unlist(GO_term))

# write all transcripts 
res_all <-
  cbind(res_all, gene_name, GO_ID, GO_term) %T>% 
  readr::write_csv(., path = file.path(outpath, "BC5_Lnc5_vs_BC6_Lnc5_diffExp_all.csv"))
###### lnc5 comparison

# get results, plot
invisible(lapply(treatments[treatments != "WT"], function(x){

  ###### get results
  # result data frame
  res_all <- 
    results(dds, contrast = c("treatment", x, "WT")) %>% 
    tibble::as_tibble(x = .) %>% 
    tibble::rownames_to_column(var = "transcript_id") %>% 
    dplyr::arrange(padj) %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::left_join(., rpkm_mean[, str_detect(colnames(rpkm_mean), paste0(x, "|WT|transcript_id"))], by = "transcript_id") %>% 
    data.table::setnames(., c(x, "WT"), c(paste0(x, "_RPKM"), "WT_RPKM")) %>% 
    dplyr::mutate(gene_symbol = mapIds(EnsDb.Mmusculus.v79, keys = transcript_id, column = "SYMBOL", keytype = "TXID", multiVals = "first"),
                  gene_biotype = mapIds(EnsDb.Mmusculus.v79, keys = transcript_id, column = "GENEBIOTYPE", keytype = "TXID", multiVals = "first"))
          
  # add gene name
  gene_name <- mapIds(org.Mm.eg.db, keys = res_all$gene_symbol, column = "GENENAME", keytype = "SYMBOL", multiVals = "first")
  gene_name[sapply(gene_name, is.null)] <- NA
  gene_name <- unname(unlist(gene_name))
  
  # add GO ID
  GO_ID <- mapIds(org.Mm.eg.db, keys = res_all$gene_symbol, column = "GO", keytype = "SYMBOL", multiVals = "first") 
  GO_ID[sapply(GO_ID, is.null)] <- NA
  GO_ID <- unname(unlist(GO_ID))
  
  # add GO term
  GO_term <- mapIds(GO.db, keys = GO_ID, column = "TERM", keytype = "GOID", multiVals = "first")
  GO_term[sapply(GO_term, is.null)] <- NA
  GO_term <- unname(unlist(GO_term))
  
  # write all transcripts 
  res_all <-
    cbind(res_all, gene_name, GO_ID, GO_term) %T>% 
    readr::write_csv(., path = file.path(outpath, str_c(x, "_vs_WT_diffExp_all.csv"))) 
  
  # write siginificantly up/downregulated transcripts
  res_sign <- 
    res_all %>% 
    dplyr::filter(padj <= 0.1) %T>% 
    readr::write_csv(., path = file.path(outpath, str_c(x, "_vs_WT_diffExp_significant.csv")))
  
  # ###### dot plot
  # res_plot <- 
  #   res_all %>% 
  #   dplyr::select(which(str_detect(colnames(.), "transcript_id|padj|log2FoldChange|RPKM"))) %>% 
  #   dplyr::mutate(sign = ifelse(padj <= 0.1, "yes", "no"),
  #                 regulation = ifelse(log2FoldChange > 0, "up", "down"), 
  #                 regulation = replace(regulation, sign == "no", NA)) %>% 
  #   dplyr::arrange(desc(padj))
  # 
  # ggplot(data = res_plot, aes(x = res_plot[, 4], y = res_plot[, 5])) +
  #   geom_point(aes(color = regulation)) +
  #   scale_x_continuous(limits = c(0, 2500)) + 
  #   scale_y_continuous(limits = c(0, 2500)) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
  #         panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank()) +
  #   ggsave(filename = file.path(outpath, str_c(x, "_vs_WT_dotplot.pdf")))
  # 
  
  ###### MA plot
  # data for plot
  res_plot <- 
    results(dds, contrast = c("treatment", x, "WT")) %>% 
    tibble::as_tibble(x = .) %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  sign = ifelse(padj < 0.1, "yes", "no"),
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, sign == "no", "not_sign"))
 # plot
 ggplot(data = res_plot, aes(x = mean, y = lfc, color = regulation, alpha = regulation)) + 
    geom_point(size = 0.1) +
    scale_x_log10(limits = c(1e-01, 1e5)) +
    scale_y_continuous(limits = c(-4, 4)) + 
    scale_colour_manual(values = c(not_sign = "gray32", down = "red3", up = "blue3")) +
    scale_alpha_manual(values = c(not_sign = 0.3, down = 1, up = 1)) +
    guides(color = FALSE, alpha = FALSE) +
    xlab("mean expression") + 
    ylab("log2FoldChange KO/WT") +
    ggsave(filename = file.path(outpath, str_c(x, "_vs_WT_MAplot.pdf")))
  
 ###### GO enrichment
 # label significantly up/downregulated transcripts 
 res_GO <- 
   results(dds, contrast = c("treatment", x, "WT")) %>% 
   as_tibble(.) %>% 
   tibble::rownames_to_column(var = "transcript_id") %>% 
   dplyr::arrange(padj) %>% 
   dplyr::filter(complete.cases(.)) %>% 
   dplyr::left_join(., trans_width, by = "transcript_id") 
 
 # get GO enrichment separately for up/down regulated genes and for all genes together
 invisible(lapply(c("upregulated", "downregulated"), function(y) {
   
   y <- "upregulated"
   
   ### set significantly up/down-regulated to 1, all else to 0 
   if(y == "upregulated"){ 
     res_regulated <- 
       res_GO %>% 
       dplyr::mutate(diff_exp = ifelse((padj <= 0.1 & log2FoldChange > 0), 1, 0))
   }else{
     if(y == "downregulated"){ 
       res_regulated <- 
         res_GO %>% 
         dplyr::mutate(diff_exp = ifelse((padj <= 0.1 & log2FoldChange < 0), 1, 0))
     }else{ # genes which are either significantly down-regulated or up-regulated
       res_regulated <- 
         res_GO %>% 
         dplyr::mutate(diff_exp = ifelse(padj <= 0.1,  1, 0))
     }
   }
   
   # probability weighting function for regulated genes
   pwf <- nullp(DEgenes = res_regulated$diff_exp %>% set_names(., res_regulated$gene_id), 
                bias.data	= res_regulated$width, 
                plot.fit = F)

   # GO enriched terms
   GO_enriched <- 
     goseq(pwf, "mm10", "ensGene") %>% 
     dplyr::mutate(p.adjust = p.adjust(over_represented_pvalue, method = "BH")) %>% 
     dplyr::select(c(1, 6, 7, 2, 3, 8, 4, 5)) %T>%
     readr::write_csv(., path = file.path(outpath, str_c(x, "_vs_WT_GOenriched_", y, "_all.csv"))) %>% 
     dplyr::filter(p.adjust < 0.05) %T>%
     readr::write_csv(., path = file.path(outpath, str_c(x, "_vs_WT_GOenriched_", y, "_significant.csv")))
     
   if(nrow(GO_enriched) > 0){
     readr::write_csv(GO_enriched, path = file.path(outpath, str_c(x, "_vs_WT_GOenriched_", y, ".csv")))
   }

 }))
 
 ###### 
  cat("done writing", x, "to", outpath, "\n")
  
}))

