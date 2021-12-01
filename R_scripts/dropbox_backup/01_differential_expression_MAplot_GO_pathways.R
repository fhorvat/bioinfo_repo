#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 24. 5. 2017.  
### AUTHOR: Filip Horvat
### NOTE:

rm(list = ls()); gc()
# options(bitmapType = 'cairo')

################################################################################### WORKING DIRECTORY
################################################################################### 
setwd("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/differential_expression_analysis")

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
library(gage)
library(pathview)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

################################################################################### PATH VARIABLES
################################################################################### 
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10"

outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Analysis/CNOT6L_2017/differential_expression_analysis"

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/Ensembl/Ensembl_GRCm38.86.20161128.gtf.gz"

sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/CNOT6L_library_size.txt"

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
  mutate(treatment = str_extract(name, "KO|WT"), 
         stage = str_extract(name, "GV|MII|1C"), 
         group = str_c(stage, "_", treatment)) %>% 
  dplyr::select(ID, name, stage, treatment, group, library_size, sample_path) %>% 
  dplyr::slice(1:18)

# read in the ENSEMBL gene table, take protein coding genes
ensembl_gtf <- read_delim(file = mm10_gtf_path, delim = "\t", comment = "#", col_names = F, col_types = cols(.default = "c")) 

# get KEGG data for mouse, keep only signaling and metabolism pathways
kegg_sigmet <- kegg.gsets("mouse")
kegg_sigmet <- kegg_sigmet$kg.sets[kegg_sigmet$sigmet.idx]

# # knownGene TxDB
# knownGene_exons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")

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
# # counting overlaps
# bamfiles <- Rsamtools::BamFileList(sample_table$sample_path, yieldSize = 2000000)
# BiocParallel::register(BiocParallel::MulticoreParam())
# se <- GenomicAlignments::summarizeOverlaps(features = knownGene_exons,
#                                            reads = bamfiles,
#                                            mode = "Union",
#                                            singleEnd = FALSE,
#                                            ignore.strand = TRUE)
# saveRDS(se, file = file.path(outpath, "se_ensembl_gtf.RDS"))
# saveRDS(se, file = file.path(outpath, "se_knownGeneTxDB_entrezID.RDS"))

# read summarisedExperiment data
se <- readRDS(file = file.path(outpath, "se_ensembl_gtf.RDS"))
# se <- readRDS(file = file.path(outpath, "se_knownGeneTxDB_entrezID.RDS"))

# FPKM normalization
fpkm_df <- 
  as.data.frame(assay(se)) %>% 
  magrittr::set_colnames(sample_table$name) %>% 
  tibble::rownames_to_column(var = "transcript_id") %>% 
  tidyr::gather(key = name, value = counts, -transcript_id) %>% 
  dplyr::left_join(., sample_table[, c("name", "library_size")], by = "name") %>% 
  dplyr::left_join(., trans_width, by = "transcript_id") %>% 
  dplyr::mutate(library_size = round(library_size / 1E6, 2), 
                width = round(width / 1E3, 2), 
                fpm = (counts / library_size), 
                fpkm = (fpm / width)) %>% 
  dplyr::select(transcript_id, name, fpkm) %>% 
  tidyr::spread(key = name, value = fpkm)

# mean FPKM in treatment 
fpkm_mean <- 
  tidyr::gather(data = fpkm_df, key = name, value = fpkm, -transcript_id) %>% 
  dplyr::left_join(., sample_table[, c("name", "group")], by = "name") %>% 
  dplyr::group_by(transcript_id, group) %>% 
  dplyr::summarise(fpkm = round(mean(fpkm), 2)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = group, value = fpkm) %>% 
  data.table::setnames(., -1, str_c("s_", colnames(.)[2:ncol(.)])) 

# set column data
colData(se) <- DataFrame(sample_table)
se$treatment %<>% factor(levels = c("WT", "KO"))
se$group %<>% factor(levels = c("GV_WT", "GV_KO", "MII_WT", "MII_KO", "1C_WT", "1C_KO"))

# experiment design: KO vs. WT for each of 3 stages
dds <- DESeqDataSet(se, design = ~group)
# dds <- dds[rowSums(counts(dds)) > 1, ]
colnames(dds) <- sample_table$name
dds <- DESeq(dds)


###### results
invisible(sapply(unique(sample_table$stage), function(stage){
  
  # compose groups 
  groups <- c(str_c(stage, "_KO"), str_c(stage, "_WT"))
  
  # get results data.frame
  results_df <- results(dds, contrast = c("group", groups))
  
  
  ####################################  
  # results with average FPKM values and GO terms (all and significantly perturbed genes) - KO vs. WT
  #################################### 
  # data frame
  res_all <- 
    results_df %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(var = "transcript_id") %>% 
    dplyr::arrange(padj) %>% 
    dplyr::filter(complete.cases(.)) %>%
    dplyr::left_join(.,   fpkm_mean[, c("transcript_id", str_c("s_", groups))], by = "transcript_id") %>% 
    data.table::setnames(., str_c("s_", groups), str_c("s_", groups, "_FPKM")) %>% 
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
    readr::write_csv(., path = file.path(outpath, str_c(stage, "_diffExp_KOvsWT_all.csv"))) 
  
  # write siginificantly up/downregulated transcripts
  res_sign <- 
    res_all %>% 
    dplyr::filter(padj <= 0.1) %T>% 
    readr::write_csv(., path = file.path(outpath, str_c(stage, "_diffExp_KOvsWT_significant.csv")))
  
  # message
  cat(stage, "writing result tables done", "\n \n \n")
  
  
  ####################################  
  # MA plot
  #################################### 
  # data for plot
  res_plot <- 
    results_df %>% 
    as.data.frame(.) %>% 
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
    ggsave(filename = file.path(outpath, str_c(stage, "_MAplot_KOvsWT.pdf")))
  
  # message
  cat(stage, "MA plot done", "\n \n \n")
  
  
  ####################################  
  # GO enrichment
  #################################### 
  # label significantly up/downregulated transcripts 
  res_GO <- 
    results_df %>% 
    as_tibble(.) %>% 
    tibble::rownames_to_column(var = "transcript_id") %>% 
    dplyr::arrange(padj) %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::left_join(., trans_width, by = "transcript_id") 
  
  # get GO enrichment separately for up/down regulated genes and for all genes together
  invisible(lapply(c("upregulated", "downregulated"), function(regulation) {
    
    ### set significantly up/down-regulated to 1, all else to 0 
    if(regulation == "upregulated"){ 
      res_regulated <- 
        res_GO %>% 
        dplyr::mutate(diff_exp = ifelse((padj <= 0.1 & log2FoldChange > 0), 1, 0))
    }else{
      if(regulation == "downregulated"){ 
        res_regulated <- 
          res_GO %>% 
          dplyr::mutate(diff_exp = ifelse((padj <= 0.1 & log2FoldChange < 0), 1, 0))
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
      dplyr::select(c(1, 6, 7, 2, 3, 8, 4, 5)) %>%
      dplyr::filter(p.adjust < 0.05) 
    
    if(nrow(GO_enriched) > 0){
      readr::write_csv(GO_enriched, path = file.path(outpath, str_c(stage, "_GOenrichment_KOvsWT_", regulation, "_significant.csv")))
    }
    
  }))
  
  # message
  cat(stage, "GO enrichment done", "\n \n \n")
  
  
  ####################################  
  # GAGE/pathview log2FC values
  #################################### 
  # results
  res_gage <- 
    results_df %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "transcript_id") %>% 
    dplyr::left_join(., trans_width, by = "transcript_id") %>% 
    dplyr::select(-width) %>% 
    mutate(entrezID = mapIds(org.Mm.eg.db, keys = gene_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"), 
           log2FoldChange = replace(log2FoldChange, is.na(log2FoldChange), 0)) %>% 
    dplyr::filter(entrezID != "", 
                  !is.na(entrezID), 
                  !str_detect(entrezID, ";"), 
                  !duplicated(entrezID))
  
  # get log2FC vector, set entrezID as names
  res_fc <- 
    res_gage$log2FoldChange %>% 
    set_names(res_gage$entrezID) 
  
  # both directions pathways (genes in pathway can be upregulated or downregulated)
  gage(exprs = res_fc, gsets = kegg_sigmet, same.dir = T) %$%
    greater %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "pid") %>% 
    dplyr::mutate(pid = str_sub(pid, 1, 8)) %>% 
    dplyr::filter(pid != "mmu04723", 
                  q.val < 0.1,
                  !is.na(q.val)) %$% 
    pid %T>% 
    sapply(., function(pid){
      pathview(gene.data = res_fc, 
               pathway.id = pid, 
               species = "mouse",
               limit = list(gene = 8, cpd = 8),
               bins = list(gene = 30, cpd = 30),
               out.suffix = str_c(stage, "_GAGE_KOvsWT_bothDirections"))
    })
  
  # message
  cat(stage, "GAGE/pathview done", "\n \n \n")
  
  
  # final message
  cat(stage, "all done!!!", "\n \n \n")
  
  
}))


