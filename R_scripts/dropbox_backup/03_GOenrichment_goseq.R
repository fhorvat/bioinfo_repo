#!/home/students/fhorvat/R/bin/Rscript
### INFO: R Script
### DATE: 18. 4. 2017.  
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

mm10_gtf_path <- "/common/WORK/fhorvat/reference/mouse/mm10/EnsemblGenes/Ensembl_GRCm38.86.20161128.gtf.gz"

sample_table_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Data/documentation/lncKO_2016_library_size.txt"

lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"

RDS_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp/results/dds_KOvsWT.RDS"

###### lnc5 comparison
RDS_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/diffExp/results/dds_design_name.RDS"
###### lnc5 comparison

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

# read dds RDS 
dds <- readRDS(file = RDS_path)

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
  gtf_trans %>%  
  unlist(.) %>% 
  as_tibble() %>% 
  dplyr::select(gene_id, transcript_id) %>% 
  unique(.) %>% 
  left_join(sum(width(gtf_trans)) %>% 
              tibble(transcript_id = names(.), width = .))

###### lnc5 comparison
# set log2FC threshold
threshold <- 1

# get results
res_GO <- 
  results(dds, contrast = c("name", "BC5_Lnc5", "BC6_Lnc5")) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(var = "transcript_id") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::mutate(diff_exp = ifelse(abs(log2FoldChange) > threshold, 1, 0))  %>% 
  dplyr::left_join(., trans_width, by = "transcript_id") 

# get GO enrichment separately for up/down regulated genes and for all genes together, write as .csv
invisible(lapply(c("upregulated", "downregulated"), function(y) {
  
  if(y == "upregulated"){ # keeps all which are not "significantly" downregulated
    res_regulated <- 
      res_GO %>% 
      dplyr::filter(!(log2FoldChange < -threshold))
  }else{
    if(y == "downregulated"){ # keeps all which are not "significantly" upregulated
      res_regulated <- 
        res_GO %>% 
        dplyr::filter(!(log2FoldChange > threshold))
    }
  }
  
  # probability weighting function for regulated genes
  pwf <- nullp(DEgenes = res_regulated$diff_exp %>% set_names(., res_regulated$gene_id), 
               bias.data	= res_regulated$width)
  
  # GO enriched terms
  GO_enriched <- 
    goseq(pwf, "mm10", "ensGene") %>% 
    dplyr::filter(p.adjust(over_represented_pvalue, method = "BH") < 0.05) %>% 
    dplyr::select(c(1, 6, 7, 2:5))
  
  if(nrow(GO_enriched) > 0){
    readr::write_csv(GO_enriched, path = file.path(outpath, stringr::str_c("BC5_Lnc5_vs_BC6_Lnc5_GOenriched_", y, ".csv")))
  }
  
}))
###### lnc5 comparison


################################################################################### 
# GO enrichment 
invisible(lapply(treatments[treatments != "WT"], function(x){
  
  ### label significantly up/downregulated transcripts for GO enrichment
  res_GO <- 
    results(dds, contrast = c("treatment", x, "WT")) %>% 
    as.data.frame(.) %>% 
    tibble::rownames_to_column(var = "transcript_id") %>% 
    dplyr::arrange(padj) %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::left_join(., trans_width, by = "transcript_id") 
  
  # get GO enrichment separately for up/down regulated genes and for all genes together, write as .csv
  invisible(lapply(c("upregulated", "downregulated", "all"), function(y) {
    
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
      }
    }else{ # genes which are either significantly down-regulated or up-regulated
      res_regulated <- 
        res_GO %>% 
        dplyr::mutate(diff_exp = ifelse(padj <= 0.1,  1, 0))
    }
    
    # probability weighting function for regulated genes
    pwf <- nullp(DEgenes = res_regulated$diff_exp %>% set_names(., res_regulated$gene_id), 
                 bias.data	= res_regulated$width)
    
    # GO enriched terms
    GO_enriched <- 
      goseq(pwf, "mm10", "ensGene") %>% 
      dplyr::filter(p.adjust(over_represented_pvalue, method = "BH") < 0.05) %>% 
      dplyr::select(c(1, 6, 7, 2:5)) 
    
    if(nrow(GO_enriched) > 0){
      readr::write_csv(GO_enriched, path = file.path(outpath, str_c(x, "_vs_WT_GOenriched_", y, ".csv")))
    }
    
  }))
  
  cat("done writing all", x, "to", outpath, "\n")
  
}))


