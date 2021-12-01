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
library(viridis)
library(dendsort)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# sort clusters
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

######################################################## PATH VARIABLES
### in and out
# set inpath
inpath <- getwd()

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
threads <- args$threads
mapped_path <- args$mapped_path
documentation_path <- args$documentation_path
features_coordinates <- args$features_coordinates
features_name <- args$features_name
genes_info_path <- args$genes_info_path
grouping_variables <- args$grouping_variables
results_groups <- args$results_groups
protein_coding_only <- args$protein_coding_only
exploratory_analysis <- args$exploratory_analysis
interactive_plots <- args$interactive_plots
counts_path <- args$counts_path
lfc_cut <- as.numeric(args$lfc_cut)
padj_cut <- as.numeric(args$padj_cut)


experiment='DicerX_mESC.DXII'
single_end=TRUE
threads=1
mapped_path='/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2018_Mar_Jun.mESC/Data/Mapped/STAR_mm10'
documentation_path='/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2018_Mar_Jun.mESC/Data/Documentation'
features_coordinates='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.gtf.gz'
features_name='ensembl.93.GRCm38.p6.20180919.UCSCseqnames'
genes_info_path='/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/ensembl.93.GRCm38.p6.20180919.UCSCseqnames.geneInfo.csv'
grouping_variables='genotype'
results_groups='RS10,RS7 RS10,RSP RSP,RS7'
protein_coding_only='no'
exploratory_analysis='yes'
interactive_plots='yes'
counts_path='./ensembl.93.GRCm38.p6.20180919.UCSCseqnames.counts.txt'
lfc_cut=as.numeric('2')
padj_cut=as.numeric('0.05')

results_groups <- str_split(results_groups, pattern = " ") %>% unlist()
grouping_variables <- str_split(grouping_variables, pattern = " ") %>% unlist()

# create and set outpath
outpath <- file.path(getwd(), str_c("results.", features_name))
dir.create(outpath)

# sample table path
sample_table_path <- list.files(documentation_path, "DicerX_mESC.2018_Jun.sampleTable.csv", full.names = T)

# reads stats path
reads_stats_path <- list.files(mapped_path, ".*\\.stats_and_tracks\\.csv", full.names = T)

# FPKM path
fpkm_path <- list.files(inpath, str_c(features_name, "\\.FPKM\\.csv$"), full.names = T)
fpkm_mean_path <- list.files(inpath, str_c(features_name, "\\.FPKM_mean\\.csv$"), full.names = T)

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.)))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats and tracks table
reads_stats <- data.table::fread(reads_stats_path)

# read FPKM tables
fpkm_tb <- readr::read_csv(fpkm_path)
fpkm_mean_tb <- readr::read_csv(fpkm_mean_path)

# read genes info
genes_info <- readr::read_csv(genes_info_path)

######################################################## MAIN CODE
### prepare tables
# get feature coordinates
features_tb <-
  counts_tb %>%
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>%
  as.data.table(.)

# get gene_id of protein coding genes
protein_genes <-
  genes_info %>%
  dplyr::filter(gene_biotype == "protein_coding") %$%
  gene_id

# counts table
se <-
  counts_tb %>%
  dplyr::select(-c(Chr:Length)) %>%
  dplyr::rename(gene_id = Geneid) %>%
  dplyr::mutate_if(is.numeric, round, digits = 0) %>%
  as.data.frame(.) %>%
  set_rownames(., .$gene_id) %>%
  dplyr::select(-gene_id) %>%
  as.matrix(.)

# join sample table with stats and tracks
sample_table[reads_stats, on = "sample_id", `:=`(library_size = genome.mapped_minus_rDNA)]

# create vector of plotly symbols in ggplot shape order
ploty_symbols <- c("square", "circle", "cross", "x", "diamond", "square-open", "circle-open", "diamond-open")


### filter samples
# prepare sample table for DESeq colData
sample_table_dds <-
  sample_table %>%
  as.data.table(.) %>%
  .[sample_id %in% str_remove_all(colnames(se), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$"), ] %>%
  .[, c("sample_id", grouping_variables), with = F] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# optionaly take only protein coding genes for analysis
if(protein_coding_only == "yes"){
  
  # filter summarizedExperiment
  se_filt <- se[rownames(se) %in% protein_genes, ]
  
}else{
  
  # don't filter
  se_filt <- se
  
}

# filter summarizedExperiment to include only chosen stage, set colData
se_filt <- SummarizedExperiment(list(counts = se_filt))
colnames(se_filt) <- str_remove_all(colnames(se_filt), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$")
se_filt <- se_filt[, colnames(se_filt)[match(rownames(sample_table_dds), colnames(se_filt))]]

# check if colnames of assay match rownames of sample table DataFrame
if(all(colnames(se_filt) == rownames(sample_table_dds))){
  
  # set sample table as colData
  colData(se_filt) <- DataFrame(sample_table_dds)
  
}else{
  
  # stop script with warrning
  stop("Columns in assay are not matching row of sample table. Please check your data annotation")
  
}

## FPM data for plots
# data for plots = log transformed counts
fpkm_log_df <-
  fpkm_tb %>%
  dplyr::select(-coordinates) %>%
  dplyr::select(gene_id, colnames(se_filt)) %>%
  dplyr::filter_at(.vars = vars(starts_with("s_")), .vars_predicate = any_vars(. > 0.5)) %>%
  dplyr::mutate_at(.vars = vars(starts_with("s_")), .funs = list(~ log2(. + 0.1))) %>%
  dplyr::filter(gene_id %in% protein_genes) %>%
  as.data.frame(.) %>%
  tibble::column_to_rownames(., var = "gene_id") %>%
  as.matrix(.)


### DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~grouped_variables)


####### EXPLORATORY ANALYSIS
### PCA plot
# data for PCA = rlog transformed counts
rlog_df <-
  vst(dds, blind = T) %>%
  assay(.)

# calculates pca
pca <-
  rlog_df %>%
  t(.) %>%
  stats::prcomp(.)

# gets percent of variance for each principal component
percentVar <- (pca$sdev)^2 / sum((pca$sdev)^2)

# makes table for ggplot
pca_tb <-
  tibble(PC1 = pca$x[, 1],
         PC2 = pca$x[, 2],
         sample_id = colnames(rlog_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))


### plot
# create bare plot object
pca_plot <- ggplot(data = pca_tb, aes(x = PC1, y = PC2, label = sample_id))

# color = first grouping variable
pca_plot <-
  pca_plot +
  geom_point(aes_string(color = grouping_variables[1], fill = grouping_variables[1]), size = 7.5, shape = 21)

# add labels, themes and save plot
pca_plot <-
  pca_plot +
  scale_color_manual(values = c("RS7" = "#619CFF", "RS10" = "#F8766D", "RSP" = "#00BA38")) +
  scale_fill_manual(values = c("RS7" = "#619CFF", "RS10" = "#F8766D", "RSP" = "#00BA38")) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(color = F, fill = F)

# save plot
ggsave(filename = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                           "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                           "png", sep = ".")),
       plot = pca_plot, width = 20, height = 10)


# add labels
pca_plot <-
  pca_plot +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", 
                   box.padding = 0.35, point.padding = 0.5, segment.color = "grey50") +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15, angle = 90)) +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5)), 
         fill = "legend")

# save labeled plot
ggsave(filename = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                           "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                           "labeled", "png", sep = ".")),
       plot = pca_plot, width = 20, height = 10)


####### DIFFERENTIAL EXPRESSION ANALYSIS
### run main DESeq2 function
# DESeq
dds_deseq <- DESeq(dds)

### shrink results
# create list of results
results_list <- purrr::map(results_groups, function(result){
  
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
      dplyr::left_join(fpkm_mean_tb %>% dplyr::select(-one_of(sample_table_dds %>%
                                                                dplyr::filter(!(grouped_variables %in% result_clean)) %$%
                                                                grouped_variables %>%
                                                                unique(.))),
                       by = "gene_id") %>%
      dplyr::mutate(comparison = str_c(result_clean[1], "_vs_", result_clean[2])) %>%
      setnames(., old = result_clean, new = str_c(result_clean, ".FPKM"))
    
  }else{
    
    # stop script with warrning
    stop(str_c("Results group ", result, " does not exist in results table. Please check you results group input!"))
    
  }
  
}) %>%
  set_names(., results_groups)


### write results
# create results workbooks
wb_all <- openxlsx::createWorkbook()
wb_significant <- openxlsx::createWorkbook()

# loop through results
invisible(purrr::map(results_groups, function(result){
  
  ## prepare results
  # get results table
  results_df <-
    results_list[[result]] %>%
    dplyr::left_join(., features_tb, by = "gene_id")
  
  # shape result
  result_clean <-
    str_split(result, pattern = ",") %>%
    unlist(.)
  
  
  ## write all results
  # add worksheet and write data
  openxlsx::addWorksheet(wb = wb_all, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]) %>% str_sub(., 1, 31))
  openxlsx::writeData(wb = wb_all, sheet = str_c(result_clean[1], "_vs_", result_clean[2]) %>% str_sub(., 1, 31), x = results_df)
  
  
  ## write only significant results, padj <= 0.1
  # filter table
  results_df_sign <-
    results_df %>%
    dplyr::filter(padj <= 0.1)
  
  # check and write
  if(nrow(results_df_sign) > 0){
    
    # add worksheet and write data
    openxlsx::addWorksheet(wb = wb_significant, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]) %>% str_sub(., 1, 31))
    openxlsx::writeData(wb = wb_significant, sheet = str_c(result_clean[1], "_vs_", result_clean[2]) %>% str_sub(., 1, 31), x = results_df_sign)
    
  }
  
}))

# save workbooks to outdir
openxlsx::saveWorkbook(wb = wb_all,
                       file = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                                       "diffExp.DESeq2", str_c(grouping_variables, collapse = "_"),
                                                       "all_results.xlsx", sep = ".")),
                       overwrite = TRUE)
openxlsx::saveWorkbook(wb = wb_significant,
                       file = file.path(outpath, str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                                       "diffExp.DESeq2", str_c(grouping_variables, collapse = "_"),
                                                       "significant_results.xlsx", sep = ".")),
                       overwrite = TRUE)


### create static MA plots
# get axis limits
results_limits <-
  results_list %>%
  dplyr::bind_rows(.) %>%
  dplyr::summarise(x_limit = baseMean %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.),
                   y_limit = log2FoldChange %>% na.omit(.) %>% abs(.) %>% max(.) %>% ceiling(.))

# loop through results
invisible(purrr::map(results_groups, function(result){
  
  ## prepare results
  # get results table
  results_df <- results_list[[result]]
  
  # shape result
  result_clean <-
    str_split(result, pattern = ",") %>%
    unlist(.)
  
  
  ## MA plot
  # data for plot
  plot_df <-
    results_df %>%
    dplyr::select(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name) %>%
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
    geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 5, shape = 20) +
    scale_x_log10(limits = c(0.1, results_limits$x_limit),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
                       breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
    scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"),
                        values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
    scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +
    guides(color = F,
           alpha = F) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(legend.title = element_blank())
  
  # save plot
  ggsave(filename = file.path(outpath,
                              str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                    "plot", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"),
                                    str_c(result_clean[1], "_vs_", result_clean[2]),
                                    "png", sep = ".")),
         plot = ma_plot, width = 20, height = 10)
  
  
  ### add gene names as labels
  # filter data
  plot_df_labels <-
    plot_df %>%
    dplyr::mutate(padj_sign = ifelse(padj <= as.numeric(padj_cut), T, F),
                  lfc_sign = ifelse((abs(lfc) >= as.numeric(lfc_cut)), T, F)) %>%
    dplyr::filter(padj_sign & lfc_sign)
  
  # annotation table
  annotations <- tibble(xpos = Inf,
                        ypos = -Inf,
                        annotateText = str_c("label cutoff: ",
                                             "p-adjusted <= ", padj_cut,
                                             ", log2FC >= ", lfc_cut))
  
  # add labels
  ma_plot_labeled <-
    ma_plot +
    geom_text(data = plot_df_labels,
              aes(x = mean, y = lfc, label = gene_name),
              check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
              colour = "black", fontface = "plain") +
    geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
              colour = "black", fontface = "italic", size = 2.5,
              hjust = 1.03, vjust = -0.5) +
    xlab("mean expression") +
    ylab(str_c("log2 fold change: ", result_clean[1], " / ", result_clean[2], "\n") %>% str_replace_all(., "_", " ")) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("gray50", "red2", "#1a75ff"))),
           alpha = F) + 
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15, angle = 90))
  
  # save plot
  ggsave(filename = file.path(outpath,
                              str_c(ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                    "plot", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"),
                                    "labeled",
                                    str_c(result_clean[1], "_vs_", result_clean[2]),
                                    "png", sep = ".")),
         plot = ma_plot_labeled, width = 20, height = 10)
  
}))

