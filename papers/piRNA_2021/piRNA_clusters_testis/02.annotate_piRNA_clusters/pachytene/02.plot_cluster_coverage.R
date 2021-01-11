### INFO: 
### DATE: Mon Aug 31 14:57:45 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/piRNA_clusters.testis/annotate_clusters/pachytene/coverage")

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

library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# clusters path
clusters_path <- file.path(inpath, "..", "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.xlsx")

# library sizes
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq/Data/Mapped/STAR_mesAur1.new"
base_path_2 <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq.reseq/Data/Mapped/STAR_mesAur1"
library_size_paths <- list.files(c(base_path, base_path_2), "library_sizes.txt", full.names = T, recursive = T)

# coverage tables path
coverage_path <- list.files(inpath, ".*MesAur1\\.1k_pachytene_clusters\\.200730\\.coverage\\.csv", full.names = T)

######################################################## READ DATA
# read clusters tb
clusters_tb <- 
  openxlsx::read.xlsx(clusters_path) %>% 
  as_tibble(., .name_repair = "unique") %>% 
  dplyr::filter(!is.na(coordinates))

# read library size table
library_size_tb <- 
  library_size_paths %>% 
  purrr::map(., readr::read_delim,  delim = "\t", col_names = c("sample_id", "library_size")) %>% 
  bind_rows(.)

# read coverage tables 
coverage_list <- purrr::map(coverage_path, readr::read_csv)

######################################################## MAIN CODE
### tidy and prepare data
# get right library sizes (19to32nt reads)
library_size_tb %<>% 
  dplyr::mutate(read_length = str_extract(sample_id, "24to31nt|21to23nt|19to32nt"), 
                read_length = replace(read_length, is.na(read_length), "full")) %>% 
  dplyr::filter(read_length == "19to32nt") %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::select(sample_id, library_size)

# get coverage in one big table, normalize to RPM, split by clusters
coverage_clusters <-
  coverage_list %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.24to31nt")) %>% 
  dplyr::left_join(., library_size_tb, by = "sample_id") %>% 
  dplyr::mutate(library_size = (library_size / 1e6),
                rpm = round((coverage / library_size), 3)) %>% 
  dplyr::select(sample_id, cluster, pos, rpm, strand) %>% 
  split(., .$cluster)

# order clusters to match table
coverage_clusters <- coverage_clusters[names(coverage_clusters[match(clusters_tb$coordinates, names(coverage_clusters))])]

# create sample table
sample_tb <- 
  tibble(sample_id = basename(coverage_path) %>% str_remove(., "\\.24to31nt\\.MesAur1.*")) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l_WT|Mov10l_KO") %>% factor(., levels = c("Mov10l_WT", "Mov10l_KO")), 
                age = str_extract(sample_id, "13dpp|21dpp") %>% factor(., levels = c("13dpp", "21dpp")), 
                sequencing = ifelse(str_detect(sample_id, "\\.reseq\\.SE$"), "reseq", "original") %>% factor(., levels = c("original", "reseq")), 
                replicate = str_extract(sample_id, "r[0-9]") %>% str_remove(., "r") %>% as.integer(.)) %>% 
  dplyr::arrange(genotype, age, replicate, sequencing)

### plot for each cluster
coverage_plot_list <- purrr::map(coverage_clusters, function(coverage_tb){
  
  # find cluster name and position in table
  cluster_name <- unique(coverage_tb$cluster)
  cluster_position <- which(clusters_tb$coordinates == cluster_name) %>% stringr::str_pad(string = ., width = 3, pad = "0")
  cluster_name <- str_replace_all(cluster_name, " ", "_")
  
  # set minus strand RPM to negative values
  plot_tb <- 
    coverage_tb %>% 
    dplyr::mutate(rpm = ifelse(strand == "+", rpm, -rpm), 
                  sample_id = str_remove(sample_id, "\\.24to31nt$"), 
                  sample_id = factor(sample_id, levels = sample_tb$sample_id)) %>% 
    dplyr::arrange(sample_id) %>% 
    dplyr::mutate(sample_id = 
                    sample_id %>% 
                    str_remove_all(., "s_testis_Mov10l_|_F[0-9]{2}-M[0-9]{1}|\\.SE") %>% 
                    str_replace_all(., "_|\\.", " ") %>% 
                    str_replace(., "reseq", "rsq")) %>% 
    dplyr::mutate(sample_id = factor(sample_id, levels = unique(.$sample_id)))
  
  # get plot limits
  plot_lim <- plot_tb$rpm %>% abs %>% max %>% ceiling
  
  # plot
  coverage_plot <-
    ggplot(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
    geom_rect() +
    geom_hline(yintercept = 0, color = "black") +
    facet_grid(rows = vars(sample_id)) + 
    coord_cartesian(ylim = c(-plot_lim, plot_lim)) +
    guides(fill = FALSE) +
    ggtitle(unique(plot_tb$cluster)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          strip.text.y = element_text(size = 7.5))
  
  # # save
  # ggsave(plot = coverage_plot,
  #        filename = file.path(outpath, "results", str_c("MesAur1.1k_pachytene_clusters.200730", cluster_position, cluster_name, "jpg", sep = ".")),
  #        width = 5,
  #        height = 1 * length(unique(coverage_tb$sample_id)), limitsize = FALSE)
  
  # return plot
  return(coverage_plot)
  
})

# save to one pdf
pdf(file = file.path(outpath, "results", "MesAur1.1k_pachytene_clusters.200730.recalculated.fused.PS.coverage.pdf"), width = 5, height = 1 * length(unique(sample_tb$sample_id)))
invisible(purrr::map(coverage_plot_list, print))
dev.off()

