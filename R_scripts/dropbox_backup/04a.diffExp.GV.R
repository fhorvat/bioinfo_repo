### INFO: Do the differential expression analysis
### DATE: Tue Dec 11 23:06:19 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Analysis/expression/joined_with_2017")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(plotly)
library(openxlsx)
library(VennDiagram)
library(geneplotter)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set input path
inpath <- getwd()

# set output path
outpath <- file.path(getwd(), "results")

# sample table
sample_path <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO/datasets/Lnc1_KO.2018_Dec/Data/Documentation/lnc1_KO.RNAseq.2017Sep_2018Dec.sampleTable.clean.csv"

# genes info
genes_info_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.91.GRCm38.p5.20180512.UCSCseqnames.geneInfo.csv"

# summarizedExperiment path
se_path <- file.path(inpath, "Lnc1_KO.2018Dec_2017Sep.GRCm38.91.reducedExons.summarizedOverlaps.RDS")

# FPKM path
fpkm_path <- file.path(inpath, "Lnc1_KO.2018Dec_2017Sep.GRCm38.91.reducedExons.FPKM.csv")

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## READ DATA
# read sample table
sample_table <- readr::read_csv(file = sample_path) 

# read additional info about genes
genes_info <- readr::read_csv(genes_info_path)

# read summarizedExperiment
se <- readRDS(file = se_path) 

# read FPKM 
fpkm_df <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
#### PREPARE DATA ####
# mutate gene info
genes_info %<>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ")

# filter sample table
sample_table_dds <- 
  sample_table %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.SE")) %>% 
  tidyr::unite(col = "group", genotype, stage, sep = "_", remove = F) %>% 
  dplyr::filter(stage == "GV") %>% 
  as.data.frame(.) %>% 
  set_rownames(., .$sample_id) 

# get gene_id of protein coding genes
protein_genes <- 
  genes_info %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# get long FPKM table
fpkm_long <- 
  fpkm_df %>% 
  dplyr::select(-c(coordinates:gene_description)) %>% 
  tidyr::gather(key = sample_id, value = fpkm, -gene_id) %>% 
  dplyr::left_join(., sample_table_dds %>% dplyr::select(sample_id, group), by = "sample_id") %>% 
  dplyr::mutate(fpkm = replace(fpkm, is.na(fpkm), 0)) 

# get mean FPKM
fpkm_mean <- 
  fpkm_long %>% 
  dplyr::group_by(gene_id, group) %>% 
  dplyr::summarise(fpkm = round(mean(fpkm), 3)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = group, value = fpkm) %>% 
  data.table::setnames(., -1, str_c(colnames(.)[2:ncol(.)], ".avgFPKM"))   


#### DESEQ2 DIFFERENTIAL EXPRESSION ####
# change colnames
colnames(se) <- str_remove(colnames(se), "\\.SE\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam")

# set stages
stages <- c("GV")

### get all results
results_all <- purrr::map(stages, function(stage){
  
  # filter sample table
  sample_table_filt <- sample_table_dds[sample_table_dds$stage == stage, ]
  
  # filter summarizedExperiment to include only protein coding genes
  se_filt <- 
    se %>% 
    .[rownames(.) %in% protein_genes, ] %>%
    .[, colnames(.)[match(rownames(sample_table_filt), colnames(.))]]
  
  # add column data to SE
  colData(se_filt) <- DataFrame(sample_table_filt)
  
  # make DESeqDataSet
  dds <- 
    DESeqDataSet(se_filt, design = ~ experiment + genotype) %>% 
    DESeq(.)
  
  # shrink results
  dds_shrink <- lfcShrink(dds, contrast = c("genotype", "Null", "WT"))
  
  # get results
  results_df <- 
    dds_shrink %>% 
    as.data.frame(.) %>% 
    as.tibble(., rownames = "gene_id") %>% 
    dplyr::arrange(padj) %>% 
    dplyr::left_join(fpkm_mean %>% dplyr::select_at(vars("gene_id", contains(stage))), by = "gene_id") %>% 
    dplyr::left_join(genes_info, by = "gene_id") %>% 
    dplyr::select_at(vars(gene_id, gene_name, baseMean, log2FoldChange, padj, 
                          matches(".*avgFPKM"), coordinates, gene_biotype, gene_description))
  
  # return results
  return(results_df)  
  
}) %>% 
  purrr::set_names(., stages)
  

### get significant results
results_significant <- purrr::map(results_all, function(results_df){
  
  # filter by p-adjusted (padj < 0.1)
  results_df %>%
    dplyr::filter(padj <= 0.1)
  
}) %>% 
  purrr::set_names(., stages)


### save results
# save all results
write.xlsx(x = results_all, 
           file = file.path(outpath, str_c("Lnc1_KO.2018Dec_2017Sep.diffExp.GRCm38.91.protein_coding.all.xlsx")), 
           asTable = rep(TRUE, length(results_all)),
           sheetName = names(results_all), 
           colWidths = "auto")

# save significant results
write.xlsx(x = results_significant, 
           file = file.path(outpath, str_c("Lnc1_KO.2018Dec_2017Sep.diffExp.GRCm38.91.protein_coding.significant.xlsx")), 
           asTable = rep(TRUE, length(results_significant)),
           sheetName = names(results_significant), 
           colWidths = "auto")


### MA plots
# save as .png
purrr::map(names(results_all), function(stage){
  
  # create plot table
  MA_df <- 
    results_all[[stage]] %>% 
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>% 
    dplyr::mutate(padj = replace(padj, is.na(padj), 1), 
                  regulation = ifelse(lfc > 0, "up", "down"), 
                  regulation = replace(regulation, padj > 0.1, "no"), 
                  regulation = factor(regulation, levels = c("no", "down", "up"))) %>% 
    dplyr::arrange(regulation)
  
  # plot 
  MA_plot <- 
    ggplot(data = MA_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation)) + 
    geom_point(size = 3, shape = 20) +
    scale_x_log10(limits = c(0.1, 1e5), 
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
    scale_y_continuous(limits = c(-4, 4), breaks = c(-4:4)) +
    scale_colour_manual(values = c(no = "gray32", down = "#1a75ff", up = "red3")) +
    scale_alpha_manual(values = c(no = 1, down = 1, up = 1)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    xlab("Mean expression") +
    ylab(str_c("log2FoldChange ", stage, " Null vs WT")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 15, vjust = - 0.2), 
          axis.title.y = element_text(size = 15, vjust = 0.3), 
          axis.text.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("Lnc1_KO.2018Dec_2017Sep.MAplot.", stage, ".Null_vs_WT.protein_coding.png")),
         plot = MA_plot, width = 15, height = 10)
  
})

# save as .html
purrr::map(names(results_all), function(stage){
  
  # create plot table
  MA_df <- 
    results_all[[stage]] %>% 
    dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, matches("\\.avgFPKM$"))) %>%
    dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                  regulation = ifelse(lfc > 0, "up", "down"),
                  regulation = replace(regulation, padj > 0.1, "no"),
                  regulation = factor(regulation, levels = c("no", "down", "up")),
                  gene_description = str_replace(gene_description, " \\[.*", ""),
                  gene_description = replace(gene_description, is.na(gene_description), "")) %>%
    dplyr::arrange(regulation)
  
  # make interactive MA plot
  if(nrow(MA_df) > 0){
    
    interactive_ma_plot <-
      plotly::plot_ly(data = MA_df,
                      x = ~mean,
                      y = ~lfc,
                      text = ~paste("</br> log2FC: ", round(lfc, 3),
                                    "</br> mean exp.: ", round(mean, 3),
                                    str_c("</br> ", "Null ", stage, " avg. FPKM:"), get(str_c("Null_", stage, ".avgFPKM")),
                                    str_c("</br> ", "WT ", stage, " avg. FPKM:"), get(str_c("WT_", stage, ".avgFPKM")),
                                    "</br>", gene_id,
                                    "</br>", gene_name,
                                    "</br>", gene_description),
                      color = ~regulation,
                      colors = c("gray32", "#1a75ff", "red3"),
                      alpha = 0.75,
                      hoverinfo = "text") %>%
      add_markers() %>%
      layout(xaxis = list(title = "mean expression", type = "log"),
             yaxis = list(title = str_c("log2FoldChange ", stage, " Null vs WT")))
    
    # save as html widget
    htmlwidgets::saveWidget(as_widget(interactive_ma_plot),
                            file = file.path(outpath, str_c("Lnc1_KO.2018Dec_2017Sep.MAplot.", stage, ".Null_vs_WT.protein_coding.html")),
                            selfcontained = T)
    
  }
  
})


