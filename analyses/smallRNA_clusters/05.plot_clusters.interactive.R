### INFO: 
### DATE: Sat Sep 08 15:03:19 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Analysis/expression/smallRNA_clusters")

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

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# cluster expression path
cluster_path <- file.path(inpath, "clusters.DicerX_embryos.WT.KO.union.rpm.3.width.50.mean_RPM.csv")

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2019_Dec/Data/Mapped/STAR_mm10/2_original_mapping"

# sample bam path
sample_bam_path <- list.files(path = mapped_path, pattern = ".*.bam", full.names = T)[1]

# annotation RDS path
annot_path <- file.path(inpath, "l.annot.RDS")

######################################################## READ DATA
# read cluster expression
cluster_df <- readr::read_csv(cluster_path)

# read annotation
l.annot <- readRDS(file = annot_path)

######################################################## MAIN CODE
# comparison name
sample_name <- "DicerX.KO.WT"

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
    as_tibble(.) 
  
  # get subjectHits
  subject_hits <- 
    x[hits$subjectHits] %>% 
    values(.) %>% 
    as.data.frame(.) %>% 
    as_tibble(.) %>% 
    tidyr::unite(gene_id, gene_biotype, gene_id, sep = ": ")
  
  # get data.frame with annotation
  annotation_df <- bind_cols(tibble(queryHits = hits$queryHits), 
                             subject_hits)
  
}) 

# reduce - left_join hits from list
hits_joined <- 
  c(list(start = tibble(queryHits = 1:length(clust_all))), 
    hits_list) %>% 
  purrr::reduce(., left_join, by = "queryHits") %>% 
  magrittr::set_colnames(., c("queryHits", names(l.annot))) %>% 
  tidyr::pivot_longer(cols = -queryHits, names_to = "category", values_to = "hit") %>% 
  dplyr::left_join(., tibble(queryHits = 1:length(clust_all), 
                             class = clust_all$class), 
                   by = "queryHits") %>% 
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
  dplyr::mutate(class = replace(class, str_detect(coordinates, "EGFP|GL3|Luc"), "siRNA_control"), 
                hit = replace(hit, class == "siRNA_control", "siRNA_control")) %>% 
  dplyr::select(coordinates, class, hit, everything()) %T>%
  readr::write_csv(., str_replace(cluster_path, "mean_RPM.csv", "mean_RPM.detailed.csv"))


## overlap with Dicer = ENSMUSG00000041415
# get Dicer1 coordinates
dicer_gr <- l.annot[["mRNA"]][str_detect(l.annot[["mRNA"]]$gene_id, "ENSMUSG00000041415")]

# get hits overlaping with Dicer
dicer_hits <- 
  findOverlaps(clust_all, dicer_gr) %>% 
  queryHits(.)


### plot clusters
# what samples to plot
x_sample <- "KO"
y_sample <- "WT"

# table for plot
plot_df <- 
  cluster_df %>% 
  dplyr::mutate(class = factor(class, levels = c("siRNA_control", "miRNA", "TE", "mRNA", "RNA_other", "other"))) %>% 
  dplyr::filter(!(1:nrow(.) %in% dicer_hits)) %>%
  dplyr::arrange(desc(class)) %>% 
  dplyr::select(x = x_sample, y = y_sample, class, coordinates, hit) %>%
  dplyr::mutate(rpm_x = x, 
                rpm_y = y,
                x = log10(x + 1), 
                y = log10(y + 1), 
                class = factor(class, levels = rev(levels(class))))

# set filter name
filter_name <- "filtered.Dicer"

### plot interactive scatterploT
# create plot
p <-
  plotly::plot_ly(data = plot_df,
                  x = ~x,
                  y = ~y,
                  text = ~paste("</br>", x_sample, "RPM:", rpm_x, 
                                "</br>", y_sample, "RPM:", rpm_y, 
                                "</br>", coordinates,
                                "</br>", class, 
                                "</br>", str_replace_all(hit, ", ", "</br> ")),
                  color = ~class,
                  colors = c("siRNA_control" = "green", "miRNA" = "cornflowerblue", "TE" = "gray20",
                             "mRNA" = "red", "RNA_other" = "gray60",
                             "other" = "gray80"),
                  alpha = 0.75,
                  hoverinfo = "text", 
                  width = 900, 
                  height = 700) %>%
  add_markers() %>%
  layout(xaxis = list(title = x_sample),
         yaxis = list(title = y_sample))

# save as html widget
htmlwidgets::saveWidget(as_widget(p),
                        file = str_replace(cluster_path, "csv", str_c(filter_name, ".html")),
                        selfcontained = T)


