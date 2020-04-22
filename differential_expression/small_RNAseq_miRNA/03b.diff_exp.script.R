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
threads <- as.numeric(args$threads)
mapped_path <- args$mapped_path
documentation_path <- args$documentation_path
features_name <- args$features_name
grouping_variables <- args$grouping_variables
results_groups <- args$results_groups
protein_coding_only <- args$protein_coding_only
exploratory_analysis <- args$exploratory_analysis
interactive_plots <- args$interactive_plots
counts_path <- args$counts_path

# create and set outpath
outpath <- file.path(getwd(), str_c("results.", features_name))
dir.create(outpath)

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# reads stats path
reads_stats_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# FPM path
fpm_path <- list.files(inpath, str_c(features_name, "\\.FPM\\.csv$"), full.names = T)
fpm_mean_path <- list.files(inpath, str_c(features_name, "\\.FPM_mean\\.csv$"), full.names = T)

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.))) %>% 
  dplyr::mutate(Geneid = make.unique(Geneid))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read stats table
reads_stats <- readr::read_delim(reads_stats_path, delim = "\t", col_names = c("sample_id", "library_size"))

# read FPM tables
fpm_tb <- readr::read_csv(fpm_path)
fpm_mean_tb <- readr::read_csv(fpm_mean_path)

######################################################## MAIN CODE
#### prepare data
## features
# get feature coordinates
features_tb <-
  counts_tb %>%
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>%
  as.data.table(.)

## sample table
# prepare sample table for DESeq colData
sample_table_dds <-
  sample_table %>%
  as.data.table(.) %>%
  .[sample_id %in% str_remove_all(colnames(counts_tb), "\\.24to31nt|\\.21to23nt|\\.genome\\.Aligned\\.sortedByCoord\\.out\\.bam|\\.total\\.bam|\\.bam$"), ] %>%
  .[, c("sample_id", grouping_variables), with = F] %>%
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>%
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

## read stats for library sizes
# get only 21-23nt reads
reads_stats %<>%
  dplyr::filter(!is.na(sample_id),
                str_detect(sample_id, "21to23nt$")) %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.21to23nt$")) %>%
  as.data.table(.)

## summarized experiment
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

# filter summarizedExperiment to include only chosen stage, set colData
se_filt <- se
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
log_df <-
  fpm_tb %>% 
  dplyr::select(-coordinates) %>% 
  dplyr::filter_at(.vars = vars(starts_with("s_")), .vars_predicate = any_vars(. > 0.5)) %>% 
  dplyr::mutate_at(.vars = vars(starts_with("s_")), .funs = list(~ log2(. + 0.1))) %>% 
  as.data.frame(.) %>% 
  tibble::column_to_rownames(., var = "gene_id") %>% 
  as.matrix(.)

## DDS
# make DESeqDataSet
dds <- DESeqDataSet(se_filt, design = ~grouped_variables)


#### exploratory analysis
### PCA plot 
## calculate
# data for PCA = rlog transformed counts
rlog_df <-
  rlog(dds, blind = T) %>%
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
         sample_id = colnames(log_df)) %>%
  dplyr::left_join(sample_table_dds , by = "sample_id") %>%
  dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " "))

## plot
# create bare plot object
pca_plot <- ggplot(data = pca_tb, aes(x = PC1, y = PC2, label = sample_id))

# if there is only one grouping variable use only color, if there is more use also a shape
if(length(grouping_variables) == 1){
  
  # color = first grouping variable
  pca_plot <-
    pca_plot +
    geom_point(aes_string(color = grouping_variables[1], fill = grouping_variables[1]), size = 5, shape = 21) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5)))
  
}else{
  
  # color = first grouping variable, shape = second grouping variable
  pca_plot <-
    pca_plot +
    geom_point(aes_string(color = grouping_variables[1], fill = grouping_variables[1], shape = grouping_variables[2]), size = 7.5) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5)),
           shape = guide_legend(override.aes = list(size = 5)))
  
}

# add labels, themes and save plot
pca_plot <-
  pca_plot +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust = - 0.2),
        axis.title.y = element_text(size = 15, vjust = - 0.2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save plot
ggsave(filename = file.path(outpath, str_c("miRBase",
                                           "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                           "png", sep = ".")),
       plot = pca_plot, width = 12, height = 10)

# add labels
pca_plot <-
  pca_plot +
  geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")

# save labeled plot
ggsave(filename = file.path(outpath, str_c("miRBase",
                                           "plot", "PCA.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                           "labeled", "png", sep = ".")),
       plot = pca_plot, width = 12, height = 10)


### distance heatmap
# calculate distance
dist_df <-
  rlog_df %>%
  t(.) %>%
  dist(.)

# make matrix
dist_matrix <- as.matrix(dist_df)
colnames(dist_matrix) <- NULL

# annotation data.frame
annotation_df <- 
  sample_table_dds %>% 
  dplyr::select(-c(sample_id, grouped_variables))

# rownames annotation
annotation_rownames <- 
  rownames(dist_matrix) %>% 
  str_remove_all(., "^s_|\\.PE$|\\.SE$") %>% 
  str_replace_all(., "_", " ")

# plot
pheatmap::pheatmap(dist_matrix,
                   clustering_distance_rows = dist_df,
                   clustering_distance_cols = dist_df,
                   col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   annotation_row = annotation_df,
                   labels_row = annotation_rownames,
                   file = file.path(outpath,
                                    str_c("miRBase",
                                          "plot", "dist_heatmap", "rlog", str_c(grouping_variables, collapse = "_"),
                                          "labeled", "png", sep = ".")),
                   height = 10,
                   width = 14)


### FPM heatmap ####
# make matrix
heatmap_matrix <- as.matrix(log_df)
rownames(heatmap_matrix) <- NULL

# annotation data.frame
annotation_df <-
  sample_table_dds %>% 
  dplyr::select(-c(sample_id, grouped_variables))

# sort rows and columns 
mat_cluster_cols <- sort_hclust(hclust(dist(t(heatmap_matrix))))
mat_cluster_rows <- sort_hclust(hclust(dist(heatmap_matrix)))

# plot
pheatmap::pheatmap(heatmap_matrix,
                   col = viridis(10),
                   annotation_col = annotation_df,
                   cluster_cols = mat_cluster_cols,
                   cluster_rows = mat_cluster_rows,
                   file = file.path(outpath,
                                    str_c("miRBase",
                                          "plot", "FPM_heatmap", "log", str_c(grouping_variables, collapse = "_"),
                                          "png", sep = ".")),
                   height = 15,
                   width = 10)



####### DIFFERENTIAL EXPRESSION ANALYSIS
# check whether to do diff. exp. analysis
if(!is.null(results_groups)){
  
  ### run main DESeq2 function
  # DESeq
  dds_deseq <- DESeq(dds)
  
  # results_groups <- str_split(results_groups, pattern = " ") %>% unlist()
  
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
        dplyr::left_join(fpm_mean_tb %>% dplyr::select(-one_of(sample_table_dds %>% 
                                                                 dplyr::filter(!(grouped_variables %in% result_clean)) %$% 
                                                                 grouped_variables %>% 
                                                                 unique(.))), 
                         by = "gene_id") %>%
        dplyr::mutate(comparison = str_c(result_clean[1], "_vs_", result_clean[2])) %>% 
        setnames(., old = result_clean, new = str_c(result_clean, ".FPM"))
      
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
    results_df <- results_list[[result]]

    # shape result
    result_clean <- 
      str_split(result, pattern = ",") %>% 
      unlist(.)
    
    
    ## write all results
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
    
  }))
  
  # save workbooks to outdir
  openxlsx::saveWorkbook(wb = wb_all, 
                         file = file.path(outpath, str_c("miRBase",
                                                         "diffExp.DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                         "all_results.xlsx", sep = ".")), 
                         overwrite = TRUE)
  openxlsx::saveWorkbook(wb = wb_significant, 
                         file = file.path(outpath, str_c("miRBase",
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
      scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit), 
                         breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
      scale_colour_manual(labels = c(no = "not significant", down = "downregulated", up = "upregulated"), 
                          values = c(no = "gray50", up = "red2", down = "#1a75ff")) +
      scale_alpha_manual(values = c(no = 0.5, down = 1, up = 1)) +      
      guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("gray50", "red2", "#1a75ff"))), 
             alpha = F) +
      xlab("mean expression") +
      ylab(str_c("log2 fold change: ", result_clean[1], " / ", result_clean[2], "\n") %>% str_replace_all(., "_", " ")) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(axis.title.x = element_text(size = 13), 
            axis.title.y = element_text(size = 13)) +
      theme(legend.title = element_blank())
    
    # turns off axis titles and legend
    ma_plot <-
      ma_plot +
      theme(legend.position = "none") +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    
    # save plot
    ggsave(filename = file.path(outpath, 
                                str_c("miRBase", 
                                      "plot", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                      str_c(result_clean[1], "_vs_", result_clean[2]), 
                                      "png", sep = ".")),
           plot = ma_plot, width = 12, height = 10)
    
    
    ### add gene IDs as labels
    # set p-adjusted and log2FC cutoff
    padj_cut <- 0.1
    lfc_cut <- 1.5
    
    # filter data
    plot_df_labels <- 
      plot_df %>% 
      dplyr::mutate(padj_sign = ifelse(padj <= padj_cut, T, F),
                    lfc_sign = ifelse((abs(lfc) >= lfc_cut), T, F)) %>% 
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
                aes(x = mean, y = lfc, label = gene_id),
                check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
                colour = "black", fontface = "plain") +
      geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText), 
                colour = "black", fontface = "italic", size = 2.5, 
                hjust = 1.03, vjust = -0.5)
    
    
    # save plot
    ggsave(filename = file.path(outpath, 
                                str_c("miRBase", 
                                      "plot", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                      "labeled", 
                                      str_c(result_clean[1], "_vs_", result_clean[2]), 
                                      "png", sep = ".")),
           plot = ma_plot_labeled, width = 12, height = 10)
    
  }))
  
  
  ### create static Vulcano plots
  # loop through results
  invisible(purrr::map(results_groups, function(result){
    
    ## prepare results
    # get results table
    results_df <- results_list[[result]]
    
    # shape result
    result_clean <- 
      str_split(result, pattern = ",") %>% 
      unlist(.)
    
    # set p-adjusted and log2FC cutoff
    padj_cut <- 0.1
    lfc_cut <- 1.5
    
    
    ## Vulcano plot - static ggplot
    # data for plot
    plot_df <-
      results_df %>%
      dplyr::select(lfc = log2FoldChange, padj, gene_id) %>%
      dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                    padj = replace(padj, padj == 0, .Machine$double.xmin),
                    padj_sign = ifelse(padj <= padj_cut, T, F),
                    lfc_sign = ifelse((abs(lfc) > lfc_cut), T, F)) %>%
      dplyr::mutate(regulation = "no",
                    regulation = replace(regulation, lfc_sign, "fold_change"),
                    regulation = replace(regulation, padj_sign, "p_value"),
                    regulation = replace(regulation, lfc_sign & padj_sign, "p_value_fold_change"),
                    regulation = factor(regulation, levels = c("no", "fold_change", "p_value", "p_value_fold_change"))) %>%
      dplyr::arrange(regulation)
    
    # annotation table
    annotations <- tibble(xpos = Inf,
                          ypos = -Inf,
                          annotateText = str_c("label cutoff: ", 
                                               "p-adjusted <= ", padj_cut, 
                                               ", log2FC >= ", lfc_cut)) 
    
    # plot
    vulcano_plot <-
      ggplot(plot_df, aes(x = lfc, y = -log10(padj), color = regulation)) +
      geom_point(size = 1.5, shape = 19, alpha = 0.5) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "longdash", colour = "black", size = 0.4) +
      geom_hline(yintercept = -log10(padj_cut), linetype = "longdash", colour = "black", size = 0.4) +
      geom_text(data = plot_df %>% dplyr::filter(padj_sign & lfc_sign),
                aes(label = plot_df %>% dplyr::filter(padj_sign & lfc_sign) %$% gene_id),
                check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
                colour = "black", fontface = "plain") +
      geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText), 
                colour = "black", fontface = "italic", size = 2.5, 
                hjust = 1.03, vjust = -0.5) +
      scale_colour_manual(values = c(no = "gray30", fold_change = "forestgreen", p_value = "royalblue", p_value_fold_change = "red2")) +
      scale_x_continuous(limits = c(-8, 8)) +
      guides(color = F, alpha = F) +
      xlab("log2FC") +
      ylab("-log10 p-adjusted value") +
      theme_bw(base_size = 24) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(legend.position = "top")
    
    # save plot
    ggsave(filename = file.path(outpath, 
                                str_c("miRBase",
                                      "plot", "vulcano", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                      str_c(result_clean[1], "_vs_", result_clean[2]),
                                      "png", sep = ".")),
           plot = vulcano_plot, width = 10, height = 10)
    
  }))
  
  
  ### create interactive plots
  if(interactive_plots == "yes"){
    
    ## MA plots
    # loop through results
    invisible(purrr::map(results_groups, function(result){
      
      ## prepare results
      # get results table
      results_df <- results_list[[result]]
      
      # shape result
      result_clean <- 
        str_split(result, pattern = ",") %>% 
        unlist(.)
      
      
      ### MA plot - interactive plot.ly
      # data for plot
      plot_df <-
        results_df %>%
        dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, coordinates, contains(".FPM"))) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      regulation = ifelse(lfc > 0, "up", "down"),
                      regulation = replace(regulation, padj > 0.1, "no"),
                      regulation = factor(regulation, levels = c("no", "down", "up"))) %>%
        dplyr::arrange(regulation)
      
      # plot
      if(nrow(plot_df) > 0){
        
        interactive_ma_plot <-
          plotly::plot_ly(data = plot_df,
                          x = ~mean,
                          y = ~lfc,
                          text = ~paste("</br> log2FC: ", round(lfc, 3),
                                        "</br> mean exp.: ", round(mean, 3),
                                        str_c("</br> ", result_clean[1], " FPM:"), get(str_c(result_clean[1], ".FPM")),
                                        str_c("</br> ", result_clean[2], " FPM:"), get(str_c(result_clean[2], ".FPM")),
                                        "</br>", gene_id,
                                        "</br>", coordinates),
                          color = ~regulation,
                          colors = c("gray32", "#1a75ff", "red3"),
                          alpha = 0.75,
                          hoverinfo = "text") %>%
          add_markers() %>%
          layout(xaxis = list(title = "mean expression", type = "log"),
                 yaxis = list(title = "log2FoldChange"))
        
        # save as html widget
        htmlwidgets::saveWidget(plotly::as_widget(interactive_ma_plot),
                                file = file.path(outpath, 
                                                 str_c("miRBase", 
                                                       "interactive", "MA", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                       str_c(result_clean[1], "_vs_", result_clean[2]),
                                                       "html", sep = ".")), 
                                selfcontained = T)
        
      }
      
    }))
    
    
    ## vulcano plots
    # loop through results
    invisible(purrr::map(results_groups, function(result){
      
      ## prepare results
      # get results table
      results_df <- results_list[[result]]
      
      # shape result
      result_clean <- 
        str_split(result, pattern = ",") %>% 
        unlist(.)
      
      # set p-adjusted and log2FC cutoff
      padj_cut <- 0.1
      lfc_cut <- 1.5
      
      ### Vulcano plot - interactive plot_ly
      # data for plot
      plot_df <-
        results_df %>%
        dplyr::select_at(.vars = vars(lfc = log2FoldChange, padj, gene_id, coordinates, contains("FPM"))) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      padj_sign = ifelse(padj <= padj_cut, T, F),
                      lfc_sign = ifelse((abs(lfc) > lfc_cut), T, F)) %>%
        dplyr::mutate(regulation = "no",
                      regulation = replace(regulation, lfc_sign, "fold_change"),
                      regulation = replace(regulation, padj_sign, "p_value"),
                      regulation = replace(regulation, lfc_sign & padj_sign, "p_value_fold_change"),
                      regulation = factor(regulation, levels = c("no", "fold_change", "p_value", "p_value_fold_change"))) %>%
        dplyr::arrange(regulation)
      
      
      # make interactive Vulcano plot
      if(nrow(plot_df) > 0){
        
        interactive_vulcano_plot <-
          plotly::plot_ly(data = plot_df,
                          x = ~lfc,
                          y = ~-log10(padj),
                          text = ~paste("</br> log2FC: ", round(lfc, 3),
                                        "</br> p-adjusted: ", round(padj, 3),
                                        str_c("</br> ", result_clean[1], " FPM:"), get(str_c(result_clean[1], ".FPM")),
                                        str_c("</br> ", result_clean[2], " FPM:"), get(str_c(result_clean[2], ".FPM")),
                                        "</br>", gene_id,
                                        "</br>", coordinates),
                          color = ~regulation,
                          colors = c("gray30", "forestgreen", "royalblue", "red2"),
                          alpha = 0.75,
                          hoverinfo = "text") %>%
          add_markers() %>%
          layout(xaxis = list(title = "log2FC"),
                 yaxis = list(title = str_c("-log10 p-adjusted value")))
        
        # save as html widget
        htmlwidgets::saveWidget(plotly::as_widget(interactive_vulcano_plot),
                                file = file.path(outpath, 
                                                 str_c("miRBase",
                                                       "interactive", "vulcano", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                       str_c(result_clean[1], "_vs_", result_clean[2]),
                                                       "html", sep = ".")),
                                selfcontained = T)
        
      }
      
    }))
    
  }
  
}

