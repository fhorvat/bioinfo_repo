### INFO: 
### DATE: Sat Sep 08 15:03:19 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/ES_DcrTrans_2012")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(data.table)
library(GenomicRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library(plotly)
library(htmlwidgets)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))
wideScreen()

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# cluster expression path
cluster_path <- list.files(inpath, pattern = "clusters.*mean.csv", full.names = T)

# annotation RDS path
annot_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/smallRNA_clusters/l.annot.RDS"

######################################################## READ DATA
# read cluster expression
cluster_df_list <- purrr::map(cluster_path, readr::read_csv)

# read annotation
l.annot <- readRDS(file = annot_path)

######################################################## MAIN CODE
### get scatterplots of cluster RPM expression
for(n in 1:length(cluster_path)){
  
  # comparison name
  sample_name <- 
    cluster_path[n] %>% 
    basename(.) %>% 
    stringr::str_extract("DcrS|Dcr\\+_-")
  
  # get cluster data.frame
  cluster_df <- cluster_df_list[[n]]
  
  # get clusters from expression table
  clust_all <- 
    cluster_df %>% 
    tidyr::separate(coordinates, c("seqnames", "start", "end"), sep = " ") %>% 
    GenomicRanges::GRanges(.)
  
  # find hits in annotation
  hits_list <- lapply(l.annot, function(x){
    
    # get hits
    hits <- 
      suppressWarnings(findOverlaps(clust_all, x)) %>%
      as.data.frame(.) %>% 
      as.tibble(.) 
    
    # get subjectHits
    subject_hits <- 
      x[hits$subjectHits] %>% 
      values(.) %>% 
      as.data.frame(.) %>% 
      as.tibble(.) %>% 
      tidyr::unite(gene_id, gene_biotype, gene_id, sep = ": ")
    
    # get data.frame with annotation
    annotation_df <- bind_cols(tibble(queryHits = hits$queryHits), 
                               subject_hits)
    
  }) 
  
  # reduce - left_join hits from list
  hits_joined <- 
    c(list(start = tibble(queryHits = 1:length(clust_all))), hits_list) %>% 
    purrr::reduce(., left_join, by = "queryHits") %>% 
    magrittr::set_colnames(., c("queryHits", names(l.annot))) %>% 
    tidyr::gather(category, hit, -queryHits) %>%  
    dplyr::left_join(tibble(queryHits = 1:length(clust_all), 
                            class = clust_all$class), by = "queryHits") %>% 
    dplyr::filter(category == class) %>% 
    dplyr::select(-class) %>% 
    dplyr::arrange(queryHits) %>% 
    dplyr::group_by(queryHits) %>% 
    dplyr::summarise(category = str_c(unique(category), collapse = ", "), 
                     hit = str_c(unique(hit), collapse = ", ")) %>% 
    dplyr::select(hit)
  
  # add hits to data.frame
  cluster_df <- 
    bind_cols(cluster_df, hits_joined) %>% 
    dplyr::select(coordinates, class, hit, everything()) %T>%
    readr::write_csv(., file.path(outpath, str_c("clusters.ESC", "DcrO", sample_name, "union", "rpm", 3, "width", 50, "mean.detailed", "csv", sep = ".")))
  
  ## overlap with Dicer = ENSMUSG00000041415
  # get Dicer1 coordinates
  dicer_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000041415")]
  
  # get hits overlaping with Dicer
  dicer_hits <- 
    findOverlaps(clust_all, dicer_gr) %>% 
    queryHits(.)
  
  ## overlap with Mos = ENSMUSG00000078365
  # get Mos genomic coordinates
  mos_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000078365")]
  
  # get hits overlaping with genomic Mos
  mos_hits <-
    findOverlaps(clust_all, mos_gr) %>%
    queryHits(.)
  
  
  ### plot clusters
  # what samples to plot
  x <- "Dcr_Rsc"
  if(sample_name == "DcrS"){
    y <- "DcrSom_Rsc"
  }else{
    y <- "Dcr+_-.tran"
  }
  
  # prepare vectors to plot (DcrO.tran, DcrS.tran, WT.utran)
  mean.x <- str_c(x, ".rpm_mean")
  mean.y <- str_c(y, ".rpm_mean")
  
  # table for plot
  plot_df <- 
    cluster_df %>% 
    dplyr::mutate(class = ifelse((1:nrow(.) %in% mos_hits), "Mos", class)) %>% 
    dplyr::mutate(class = factor(class, levels = c("Mos", "miRNA", "TE", "mRNA", "RNA_other", "other"))) %>% 
    dplyr::filter(!(1:nrow(.) %in% dicer_hits)) %>%
    dplyr::filter(!(str_detect(coordinates, "pCAG-EGFP_MosIR"))) %>%
    dplyr::arrange(desc(class)) %>% 
    dplyr::select(x = mean.x, y = mean.y, class, coordinates, hit) %>%
    dplyr::mutate(x = ifelse(class == "Mos", x * 2, x))
  
  # multiply by 2 Mos in DcrS mean sample
  if(sample_name == "DcrS"){
    plot_df %<>%
      dplyr::mutate(y = ifelse(class == "Mos", y * 2, y))
  }
  
  # log10 values for plot
  plot_df %<>% 
    dplyr::mutate(rpm_x = x, 
                  rpm_y = y,
                  x = log10(x + 1), 
                  y = log10(y + 1), 
                  class = factor(class, levels = rev(levels(class))))
  
  # set filter name
  filter_name <- "filtered.MosIR.Dicer"
  
  ### plot interactive scatterploT
  # create plot
  p <-
    plotly::plot_ly(data = plot_df,
                    x = ~x,
                    y = ~y,
                    text = ~paste("</br>", str_remove(mean.x, ".tran.*|.utran.*|.rpm_mean.*"), "RPM:", rpm_x, 
                                  "</br>", str_remove(mean.y, ".tran.*|.utran.*|.rpm_mean.*"), "RPM:", rpm_y, 
                                  "</br>", coordinates,
                                  "</br>", class, 
                                  "</br>", str_replace_all(hit, ", ", "</br> ")),
                    color = ~class,
                    colors = c("Mos" = "red", "miRNA" = "cornflowerblue", "TE" = "gray20",
                               "mRNA" = "red", "RNA_other" = "gray45",
                               "other" = "gray80"),
                    alpha = 0.75,
                    hoverinfo = "text", 
                    width = 900, 
                    height = 700) %>%
    add_markers() %>%
    layout(xaxis = list(title = x),
           yaxis = list(title = y))
  
  # save as html widget
  htmlwidgets::saveWidget(as_widget(p),
                          file = file.path(outpath, str_c("clusters.ESC", "DcrO", sample_name, "union", "rpm", 3,
                                                          "width", 50, "meanRPM", filter_name, "html", sep = ".")),
                          selfcontained = T)
  
}


