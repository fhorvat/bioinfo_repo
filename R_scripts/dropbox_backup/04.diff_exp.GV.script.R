### INFO: 
### DATE: Thu Apr 25 16:58:31 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(BiocParallel)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)
library(plotly)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

# set outpath
outpath <- file.path(getwd(), "results")

# # get arguments from command line, transform to named vector
# args <-
#   commandArgs(trailingOnly = TRUE) %>%
#   parseCommandLineArguments(.)

args <- 
  "--experiment Wu_2019_unpub_GSE133748 \
--ensembl_version 93 \
--genome_path /common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2 \
--mapped_path /common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Wu_2019_unpub_GSE133748/Data/Mapped/STAR_mm10 \
--threads 1 \
--grouping_variables stage fake_genotype \
--results_groups GV_fake_KO,GV_WT \
--protein_coding no \
--exploratory_analysis no" %>% 
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
ensembl_version <- args$ensembl_version
genome_path <- args$genome_path
mapped_path <- args$mapped_path
threads <- args$threads
grouping_variables <- args$grouping_variables
results_groups <- args$results_groups
protein_coding <- args$protein_coding
exploratory_analysis <- args$exploratory_analysis



### other experiment paths
# set base experiment path
base_path <- file.path("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq", experiment)

# documentation path
documentation_path <- file.path(base_path, "Data/Documentation")

# analysis path
analysis_path <- file.path(base_path, "Analysis/expression")


### documentation
# sample table path
sample_table_path <- list.files(path = documentation_path, pattern = ".*\\.sampleTable\\.csv", full.names = T)

# stats and tracks path
stats_and_tracks_path <- list.files(path = mapped_path, pattern = ".*\\.stats_and_tracks\\.csv", full.names = T)

# summarizedExperiment path
se_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.se\\.RDS$"), full.names = T)

# FPKM table path
fpkm_path <- list.files(path = analysis_path, pattern = str_c("ensembl\\.", ensembl_version, ".*UCSCseqnames.*\\.FPKM_mean\\.csv$"), full.names = T)


### genome
# gene info path
genes_info_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# reduced exons path
exons_path <- list.files(path = genome_path, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.reducedExons.RDS$"), full.names = T)

######################################################## READ DATA
# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats and tracks table
stats_and_tracks <- data.table::fread(stats_and_tracks_path)

# summarizedExperiment 
se <- readRDS(se_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

# read ENSEMBL reduced exons
exons_gr <- readRDS(file = exons_path)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
### prepare data
# join sample table with stats and tracks
sample_table[stats_and_tracks, on = "sample_id", `:=`(library_size = genome.mapped_minus_rDNA)]

# filter ensembl genes info
genes_info_tidy <- 
  genes_info %>% 
  dplyr::mutate(cooridnates = str_c(seqnames, ":", start, "-", end, "|", strand)) %>% 
  dplyr::select(-c(seqnames:strand))

# get gene_id of protein coding genes
protein_genes <- 
  genes_info_tidy %>% 
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# create vector of plotly symbols in ggplot shape order
ploty_symbols <- c("square", "circle", "cross", "x", "diamond", "square-open", "circle-open", "diamond-open")


### filter samples
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  dplyr::filter(stage == "GV", genotype == "WT") %>% 
  dplyr::mutate(row_n = 1:n(), 
                fake_genotype = ifelse((row_n %% 2) == 0, "WT", "fake_KO")) %>% 
  as.data.table(.) %>% 
  .[, c("sample_id", grouping_variables), with = F] %>% 
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# optionaly take only protein coding genes for analysis
if(protein_coding == "yes"){
  
  # filter summarizedExperiment
  se_filt <- se[rownames(se) %in% protein_genes, ]
  
}else{
  
  # don't filter
  se_filt <- se
  
}

# filter summarizedExperiment to include only chosen stage, set colData
colnames(se_filt) <- str_remove(colnames(se_filt), "\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$")
se_filt <- se_filt[, colnames(se_filt)[match(rownames(sample_table_dds), colnames(se_filt))]]

# check if colnames of assay match rownames of sample table DataFrame
if(all(colnames(se_filt) == rownames(sample_table_dds))){
  
  # set sample table as colData
  colData(se_filt) <- DataFrame(sample_table_dds)
  
}else{
  
  # stop script with warrning
  stop("Columns in assay are not matching row of sample table. Please check your data annotation")
  
}


### DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~grouped_variables)


####### DIFFERENTIAL EXPRESSION ANALYSIS
### run main DESeq2 function
# DESeq
dds_deseq <- DESeq(dds)


### shrink results
# create list of results
result <- results_groups[1]

# shape result
result_clean <- 
  str_split(result, pattern = ",") %>% 
  unlist(.)

# check if results groups make sense
if((length(result_clean) == 2) & (all(result_clean %in% sample_table_dds$grouped_variables))){
  
  # get results, shrink logFC
  dds_shrink <- lfcShrink(dds_deseq, contrast = c("grouped_variables", result_clean[1], result_clean[2]))
  
  # get results table, add to sheet
  results_df <-
    dds_shrink %>%
    as_tibble(., rownames = "gene_id") %>%
    dplyr::arrange(padj) %>%
    dplyr::left_join(fpkm_tb %>% dplyr::select(gene_id, gene_name:gene_description), 
                     by = "gene_id") %>%
    dplyr::mutate(comparison = str_c(result_clean[1], "_vs_", result_clean[2]))
  
}else{
  
  # stop script with warrning 
  stop(str_c("Results group ", result, " does not exist in results table. Please check you results group input!"))      
  
}



### write results
# create results workbooks
wb_all <- openxlsx::createWorkbook()
wb_significant <- openxlsx::createWorkbook()

# add worksheet and write data
openxlsx::addWorksheet(wb = wb_all, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]))
openxlsx::writeData(wb = wb_all, sheet = str_c(result_clean[1], "_vs_", result_clean[2]), x = results_df)

## write only significant results, padj <= 0.1
# filter table
results_df_sign <- 
  results_df %>% 
  dplyr::filter(padj <= 0.1)

# check and write
if(nrow(results_df_sign) > 0){
  
  # add worksheet and write data
  openxlsx::addWorksheet(wb = wb_significant, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]))
  openxlsx::writeData(wb = wb_significant, sheet = str_c(result_clean[1], "_vs_", result_clean[2]), x = results_df_sign)
  
}

# save workbooks to outdir
openxlsx::saveWorkbook(wb = wb_all, 
                       file = file.path(outpath, str_c("diffExp.DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                       ifelse(protein_coding == "yes", "protein_coding", "all_biotype"),
                                                       "all_results.xlsx", sep = ".")), 
                       overwrite = TRUE)

openxlsx::saveWorkbook(wb = wb_significant, 
                       file = file.path(outpath, str_c("diffExp.DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                       ifelse(protein_coding == "yes", "protein_coding", "all_biotype"),
                                                       "significant_results.xlsx", sep = ".")), 
                       overwrite = TRUE)


### create static MA plots
# get axis limits
results_limits <- 
  results_df %>% 
  dplyr::summarise(x_limit = baseMean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.), 
                   y_limit = log2FoldChange %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# data for plot
plot_df <-
  results_df %>%
  dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id) %>%
  dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                padj = replace(padj, padj == 0, .Machine$double.xmin),
                sign = ifelse(padj < 0.1, "yes", "no"),
                regulation = ifelse(lfc > 0, "up", "down"),
                regulation = replace(regulation, padj > 0.1, "no"),
                regulation = factor(regulation, levels = c("no", "up", "down"))) %>%
  dplyr::arrange(regulation)

# plot
ma_plot <-
  ggplot() + 
  geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 2.5, shape = 20) +
  scale_x_log10(limits = c(0.01, results_limits$x_limit),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(limits = c(-5, 5), breaks = c(-5:5)) +
  scale_colour_manual(values = c(no = "gray50", down = "#1a75ff", up = "red2")) +
  scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +      
  guides(color = FALSE, alpha = FALSE) +
  xlab("mean expression") +
  ylab("log2FC") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("diffExp.DESeq2", 
                                           str_c(grouping_variables, collapse = "_"), 
                                           str_c(result_clean[1], "_vs_", result_clean[2]),
                                           ifelse(protein_coding == "yes", "protein_coding", "all_biotype"),
                                           "MA_plot.png", sep = ".")),
       plot = ma_plot, width = 10, height = 10)



# ### create interactive MA plots
# # data for plot
# plot_df <-
#   results_df %>%
#   dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description)) %>%
#   dplyr::mutate(padj = replace(padj, is.na(padj), 1),
#                 padj = replace(padj, padj == 0, .Machine$double.xmin),
#                 regulation = ifelse(lfc > 0, "up", "down"),
#                 regulation = replace(regulation, padj > 0.1, "no"),
#                 regulation = factor(regulation, levels = c("no", "down", "up")),
#                 gene_description = str_remove(gene_description, " \\[.*"),
#                 gene_description = replace(gene_description, is.na(gene_description), "")) %>%
#   dplyr::arrange(regulation)
# 
# # plot
# if(nrow(plot_df) > 0){
#   
#   interactive_ma_plot <-
#     plotly::plot_ly(data = plot_df,
#                     x = ~mean,
#                     y = ~lfc,
#                     text = ~paste("</br> log2FC: ", round(lfc, 3),
#                                   "</br> mean exp.: ", round(mean, 3),
#                                   "</br>", gene_id,
#                                   "</br>", gene_name,
#                                   "</br>", gene_description),
#                     color = ~regulation,
#                     colors = c("gray32", "#1a75ff", "red3"),
#                     alpha = 0.75,
#                     hoverinfo = "text") %>%
#     add_markers() %>%
#     layout(xaxis = list(title = "mean expression", type = "log"),
#            yaxis = list(title = "log2FoldChange"))
#   
#   # save as html widget
#   htmlwidgets::saveWidget(plotly::as_widget(interactive_ma_plot),
#                           file = file.path(outpath, str_c("diffExp.DESeq2", 
#                                                           str_c(grouping_variables, collapse = "_"), 
#                                                           str_c(result_clean[1], "_vs_", result_clean[2]),
#                                                           ifelse(protein_coding == "yes", "protein_coding", "all_biotype"),
#                                                           "MA_plot.html", sep = ".")), 
#                           selfcontained = T)
#   
# }
