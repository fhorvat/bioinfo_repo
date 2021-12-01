### INFO: 
### DATE: Thu Jul 18 19:52:56 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/miRNA.Shubha/Analysis/2020_paper/Freimer_microarray/highlight_mir15a_targets")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
### in and out
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()


## gene annotation
# set ensembl version
ensembl_version <- 93

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# gene info path
genes_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.geneInfo.csv$"), full.names = T)

# transcripts info path
transcripts_info_path <- list.files(path = genome_dir, pattern = str_c("ensembl.", ensembl_version, ".*UCSCseqnames.transcriptInfo.csv$"), full.names = T)


## results
# diff. expression results path
results_path <- file.path(inpath, "..", "Freimer_2017.diff_exp.CrePos_vs_CreNeg.miR-15a.csv")

# targets list path
targets_list_path <- list.files(path = inpath, pattern = ".*miR_15a_targets\\.csv", full.names = T)

######################################################## READ DATA
# read results
results_df <- readr::read_csv(results_path)

# significant results
results_df_sign <-
  results_df %>%
  dplyr::filter(as.numeric(padj) <= 0.05) %>% 
  dplyr::mutate(regulation = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::select(ID, regulation)

# read targets
targets_list <- 
  purrr::map(targets_list_path, readr::read_csv) %>% 
  set_names(., str_extract(targets_list_path, "miRBase|miRDB"))

######################################################## MAIN CODE
### plot MA plots
purrr::map(names(targets_list), function(target_base){
  
  # get targets
  targets_tb <- targets_list[[target_base]]
  
  # get top 200 targets
  targets_top <- 
    targets_tb %>% 
    dplyr::inner_join(., results_df, by = "gene_name") %>% 
    dplyr::select(ID, rank) %>% 
    dplyr::mutate(counter = group_indices(., rank)) %>% 
    dplyr::filter(counter <= 200) %>% 
    dplyr::mutate(regulation = "top_targets") %>% 
    dplyr::select(ID, regulation)

  # data for plot
  plot_df <-
    results_df %>%
    dplyr::select(mean = log2_mean_exp, lfc = log2FoldChange, padj, ID, gene_name) %>%
    dplyr::left_join(., targets_top, by = "ID") %>%
    dplyr::mutate(regulation = replace(regulation, is.na(regulation), "other"),
                  regulation = factor(regulation, levels = c("other", "top_targets"))) %>%
    dplyr::arrange(regulation)
  
  # labels data
  plot_df_labels <-
    plot_df %>%
    dplyr::filter(regulation == "top_targets") %>%
    dplyr::mutate(gene_name = ifelse(is.na(gene_name), ID, gene_name))
  
  # annotation table
  annotations <- tibble(xpos = Inf,
                        ypos = -Inf,
                        annotateText = str_c("label cutoff: ", target_base,
                                             " top 200 targets"))
  
  # plot
  ma_plot <-
    ggplot() +
    geom_point(data = plot_df, aes(x = mean, y = lfc, color = regulation, alpha = regulation), size = 5, shape = 20) +
    # scale_x_log10(limits = c(0.01, results_limits$x_limit),
    #               breaks = scales::trans_breaks("log10", function(x) 10^x),
    #               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # scale_y_continuous(limits = c(-results_limits$y_limit, results_limits$y_limit),
    #                    breaks = c(-results_limits$y_limit:results_limits$y_limit)) +
    scale_colour_manual(labels = c(other = "other", top_targets = "top targets"),
                        values = c(other = "gray50", top_targets = "red2")) +
    scale_alpha_manual(values = c(other = 0.5, top_targets = 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  # add labels
  ma_plot_labeled <-
    ma_plot +
    # geom_text(data = plot_df_labels,
    #           aes(x = mean, y = lfc, label = gene_name),
    #           check_overlap = TRUE, size = 3, hjust = 0, vjust = 1.5,
    #           colour = "black", fontface = "plain") +
    # geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText),
    #           colour = "black", fontface = "italic", size = 2.5,
    #           hjust = 1.03, vjust = -0.5) +
    guides(color = guide_legend(override.aes = list(shape = 23, size = 5, fill = c("gray50", "red2"))),
           alpha = F) +
    xlab("log2 mean expression") +
    ylab(str_c("log2 fold change: ", "Cre pos.", " / ", "Cre neg.", " miR-15a", "\n") %>% str_replace_all(., "_", " ")) +
    theme(axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13, angle = 90)) +
    theme(legend.position = "bottom")
  
  # save plot
  ggsave(filename = file.path(outpath, str_c("Freimer_2017", "MA_plot", "CrePos_vs_CreNeg", "miR-15a", 
                                             str_c(target_base, "_top200_targets"), "png", sep = ".")),
         plot = ma_plot_labeled, width = 10, height = 10)

})
