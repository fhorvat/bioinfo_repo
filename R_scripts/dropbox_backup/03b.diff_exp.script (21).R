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

# get arguments from command line, transform to named vector
args <-
  commandArgs(trailingOnly = TRUE) %>%
  parseCommandLineArguments(.)

# arguments from command line
experiment <- args$experiment
ensembl_version <- args$ensembl_version
genome_path <- args$genome_path
mapped_path <- args$mapped_path
threads <- args$threads
grouping_variables <- args$grouping_variables
results_groups <- args$results_groups
protein_coding_only <- args$protein_coding_only
exploratory_analysis <- args$exploratory_analysis
interactive <- args$interactive

experiment='hamster_oocyte_Mov10l.RNAseq'
single_end=TRUE
threads=1
mapped_path='/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Mapped/bbmap_mesAur1'
documentation_path='/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Data/Documentation'
feature_coordinates='/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.RNAseq/Analysis/expression_coverage/joined_coverage.perfect_multimappers.coordinates.saf'
feature_names='joined_coverage.perfect_multimappers.coordinates'
grouping_variables='genotype'
results_groups='--results_groups Mov10l_KO,Mov10l_WT Mov10l_HET,Mov10l_WT Mov10l_KO,Mov10l_HET' %>% parseCommandLineArguments(.)
protein_coding_only='no'
exploratory_analysis='yes'
interactive_plots='yes'


# counts path
counts_path <- list.files(inpath, ".*\\.counts\\.txt$", full.names = T)

# sample table path
sample_table_path <- list.files(documentation_path, ".*sampleTable\\.csv", full.names = T)

# reads stats path
read_stats_path <- file.path(mapped_path, "3_logs", "log.read_stats.txt")

# FPKM path
fpkm_path <- list.files(inpath, ".*\\.FPKM_mean\\.csv$", full.names = T)

######################################################## READ DATA
# read counts from featureCounts
counts_tb <-
  readr::read_delim(counts_path, delim = "\t", comment = "#") %>%
  set_colnames(., basename(colnames(.)))

# read sample table
sample_table <- data.table::fread(sample_table_path)

# read and clean stats
reads_stats <-
  read_delim(read_stats_path, delim = "\t") %>%
  dplyr::select(sample_id, library_size = genome.mapped_minus_rDNA) %>% 
  as.data.table(.)

# read FPKM table
fpkm_tb <- readr::read_csv(fpkm_path)

######################################################## MAIN CODE
### prepare tables
# get feature coordinates
features_tb <- 
  counts_tb %>% 
  dplyr::select(gene_id = Geneid, seqnames = Chr, start = Start, end = End, width = Length) %>% 
  as.data.table(.)

# counts table
se <- 
  counts_tb %>%
  dplyr::select(-c(Chr:Length)) %>% 
  dplyr::rename(gene_id = Geneid) %>% 
  as.data.frame(.) %>% 
  set_rownames(., .$gene_id) %>% 
  dplyr::select(-gene_id) %>% 
  as.matrix(.)

# join sample table with stats and tracks
sample_table[reads_stats, on = "sample_id", `:=`(library_size = library_size)]

# create vector of plotly symbols in ggplot shape order
ploty_symbols <- c("square", "circle", "cross", "x", "diamond", "square-open", "circle-open", "diamond-open")


### filter samples
# prepare sample table for DESeq colData
sample_table_dds <- 
  sample_table %>% 
  as.data.table(.) %>% 
  .[, c("sample_id", "sample_name", grouping_variables), with = F] %>% 
  .[, grouped_variables := do.call(str_c, c(.SD, sep = "_")), .SDcols = grouping_variables] %>% 
  as.data.frame(.) %>%
  set_rownames(., .$sample_id)

# filter summarizedExperiment to include only chosen stage, set colData
se_filt <- SummarizedExperiment(list(counts = se))
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


####### EXPLORATORY ANALYSIS 
if(exploratory_analysis == "yes"){
  
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
    dplyr::mutate(sample_id = str_replace(sample_id, "s_|r", "") %>% str_replace_all(., "_", " ") %>% str_c(., " ", sample_name))
 
  
  ### plot
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
  ggsave(filename = file.path(outpath, str_c("PCA_plot.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                             ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "png", sep = ".")),
         plot = pca_plot, width = 12, height = 10)
  
  # add labels
  pca_plot <- 
    pca_plot +
    geom_label_repel(aes(label = sample_id), fontface = "bold", color = "black", box.padding = 0.35, point.padding = 0.5, segment.color = "grey50")
  
  # save labeled plot
  ggsave(filename = file.path(outpath, str_c("PCA_plot.PC1_PC2", "rlog", str_c(grouping_variables, collapse = "_"),
                                             ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "labeled", "png", sep = ".")),
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
    dplyr::select(-c(sample_id, sample_name, grouped_variables))
  
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
                                      str_c("dist_heatmap", "rlog", str_c(grouping_variables, collapse = "_"),
                                            ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "png", sep = ".")),
                     height = 10,
                     width = 14)
  
}


####### DIFFERENTIAL EXPRESSION ANALYSIS
# check whether to do diff. exp. analysis
if(!is.null(results_groups)){
  
  ### run main DESeq2 function
  # DESeq
  dds_deseq <- DESeq(dds)
  
  ### shrink results
  # create list of results
  results_list <- purrr::map(results_groups[[1]], function(result){
    
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
        dplyr::left_join(fpkm_tb %>% dplyr::select(-one_of(sample_table_dds %>% dplyr::filter(!(grouped_variables %in% result_clean)) %$% grouped_variables %>% unique(.))), 
                         by = "gene_id") %>%
        dplyr::mutate(comparison = str_c(result_clean[1], "_vs_", result_clean[2])) %>% 
        setnames(., old = result_clean, new = str_c(result_clean, ".FPKM"))
      
    }else{
      
      # stop script with warrning 
      stop(str_c("Results group ", result, " does not exist in results table. Please check you results group input!"))      
      
    }
    
  }) %>% 
    set_names(., results_groups[[1]])
  
  
  ### write results
  # create results workbooks
  wb_all <- openxlsx::createWorkbook()
  wb_significant <- openxlsx::createWorkbook()
  
  # loop through results
  invisible(purrr::map(results_groups[[1]], function(result){
    
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
      dplyr::filter(padj <= 0.1) %>% 
      dplyr::mutate(coordinates = gene_id) %>% 
      tidyr::separate(coordinates, into = c("seqnames", "coordinates"), sep = ":") %>% 
      tidyr::separate(coordinates, into = c("start", "end"), sep = "-") %>% 
      dplyr::mutate(width = as.numeric(end) - as.numeric(start))
    
    # check and write
    if(nrow(results_df_sign) > 0){
      
      # add worksheet and write data
      openxlsx::addWorksheet(wb = wb_significant, sheetName = str_c(result_clean[1], "_vs_", result_clean[2]))
      openxlsx::writeData(wb = wb_significant, sheet = str_c(result_clean[1], "_vs_", result_clean[2]), x = results_df_sign)
      
    }
    
  }))
  
  # save workbooks to outdir
  openxlsx::saveWorkbook(wb = wb_all, 
                         file = file.path(outpath, str_c("diffExp.DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                         ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
                                                         "all_results.xlsx", sep = ".")), 
                         overwrite = TRUE)
  openxlsx::saveWorkbook(wb = wb_significant, 
                         file = file.path(outpath, str_c("diffExp.DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                         ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"),
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
  invisible(purrr::map(results_groups[[1]], function(result){
    
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
                                str_c("MA_plot", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                      str_c(result_clean[1], "_vs_", result_clean[2]),
                                      ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "png", 
                                      sep = ".")),
           plot = ma_plot, width = 12, height = 10)
    
  }))
  
  
  ### create interactive MA plots
  if(interactive == "yes"){
    
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
        dplyr::select_at(.vars = vars(mean = baseMean, lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, contains("FPKM"))) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      regulation = ifelse(lfc > 0, "up", "down"),
                      regulation = replace(regulation, padj > 0.1, "no"),
                      regulation = factor(regulation, levels = c("no", "down", "up")),
                      gene_description = str_remove(gene_description, " \\[.*"),
                      gene_description = replace(gene_description, is.na(gene_description), "")) %>%
        dplyr::arrange(regulation)
      
      # plot
      if(nrow(plot_df) > 0){
        
        interactive_ma_plot <-
          plotly::plot_ly(data = plot_df,
                          x = ~mean,
                          y = ~lfc,
                          text = ~paste("</br> log2FC: ", round(lfc, 3),
                                        "</br> mean exp.: ", round(mean, 3),
                                        str_c("</br> ", result_clean[1], " FPKM:"), get(str_c(result_clean[1], ".FPKM")),
                                        str_c("</br> ", result_clean[2], " FPKM:"), get(str_c(result_clean[2], ".FPKM")),
                                        "</br>", gene_id,
                                        "</br>", gene_name,
                                        "</br>", gene_description),
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
                                                 str_c("MA_plot", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                       str_c(result_clean[1], "_vs_", result_clean[2]),
                                                       ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "html", 
                                                       sep = ".")), 
                                selfcontained = T)
        
      }
      
    }))
    
  }
  
  
  
  ### create static Vulcano plots
  # loop through results
  invisible(purrr::map(results_groups[[1]], function(result){
    
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
                                str_c("vulcano_plot", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                      str_c(result_clean[1], "_vs_", result_clean[2]),
                                      ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "png", 
                                      sep = ".")),
           plot = vulcano_plot, width = 10, height = 10)
    
  }))
  
  
  
  ### create interactive Vulcano plots
  if(interactive == "yes"){
    
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
        dplyr::select_at(.vars = vars(lfc = log2FoldChange, padj, gene_id, gene_name, gene_description, contains("FPKM"))) %>%
        dplyr::mutate(padj = replace(padj, is.na(padj), 1),
                      padj = replace(padj, padj == 0, .Machine$double.xmin),
                      padj_sign = ifelse(padj <= padj_cut, T, F),
                      lfc_sign = ifelse((abs(lfc) > lfc_cut), T, F)) %>%
        dplyr::mutate(regulation = "no",
                      regulation = replace(regulation, lfc_sign, "fold_change"),
                      regulation = replace(regulation, padj_sign, "p_value"),
                      regulation = replace(regulation, lfc_sign & padj_sign, "p_value_fold_change"),
                      regulation = factor(regulation, levels = c("no", "fold_change", "p_value", "p_value_fold_change"))) %>%
        dplyr::mutate(gene_description = str_remove(gene_description, " \\[.*"),
                      gene_description = replace(gene_description, is.na(gene_description), "")) %>%
        dplyr::arrange(regulation)
      
      
      # make interactive Vulcano plot
      if(nrow(plot_df) > 0){
        
        interactive_vulcano_plot <-
          plotly::plot_ly(data = plot_df,
                          x = ~lfc,
                          y = ~-log10(padj),
                          text = ~paste("</br> log2FC: ", round(lfc, 3),
                                        "</br> p-adjusted: ", round(padj, 3),
                                        str_c("</br> ", result_clean[1], " FPKM:"), get(str_c(result_clean[1], ".FPKM")),
                                        str_c("</br> ", result_clean[2], " FPKM:"), get(str_c(result_clean[2], ".FPKM")),
                                        "</br>", gene_id,
                                        "</br>", gene_name,
                                        "</br>", gene_description),
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
                                                 str_c("vulcano_plot", "DESeq2", str_c(grouping_variables, collapse = "_"), 
                                                       str_c(result_clean[1], "_vs_", result_clean[2]),
                                                       ifelse(protein_coding_only == "yes", "protein_coding", "all_biotype"), "html", 
                                                       sep = ".")),
                                selfcontained = T)
        
      }
      
    }))
    
  }
  
}
