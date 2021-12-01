#!/common/WORK/fhorvat/programi/R/R/bin/Rscript
### INFO: combine track URLs and read stats table
### DATE: 09. 10. 2017.
### AUTHOR: Filip Horvat

rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/RNAseq/Ye_2017_RNA_GSE95145/Data/Mapped/STAR_mm10")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(ggplot2)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
outpath <- getwd()
stats_path <- file.path(outpath, "log.read_stats.txt")
tracks_path <- file.path(outpath, "log.tracks_URL.txt")

experiment <- "Ye_2017_RNA_GSE95145"
experiment_name <- "Ye_2017"

######################################################## READ DATA

######################################################## MAIN CODE
# stats
stats_tbl <- readr::read_delim(stats_path, delim = "\t")

# tracks
tracks_tbl <-
  readr::read_delim(file = tracks_path, col_names = F, delim = "\t") %>%
  dplyr::select(URL = X1) %>%
  dplyr::filter(str_detect(string = URL, pattern = "\\.bw$|\\.bam$"),
                str_detect(string = URL, pattern = "http"),
                !str_detect(string = URL, pattern = "bai"),
                str_detect(string = URL, pattern = str_c(stats_tbl$sample_id, collapse = "|"))) %>%
  dplyr::mutate(sample_id = basename(URL) %>% str_remove_all(., "\\.bam$|\\.bw$|\\.scaled"),
                experiment = experiment,
                scaled = ifelse(str_detect(URL, "scaled"), "RPM_scaled", "raw"),
                file_type = ifelse(test = str_detect(basename(URL), "bw"),
                                   yes = "coverage",
                                   no = "individual_reads"),
                bw_name = ifelse(test = (scaled == "RPM_scaled"),
                                 yes = str_c(str_remove(sample_id, "^s_"), experiment_name, "scaled", sep = "."),
                                 no = str_c(str_remove(sample_id, "^s_"), experiment_name, "raw", sep = ".")),
                bam_name = str_c(str_remove(sample_id, "^s_"), experiment_name, "bam", sep = "."),
                URL = ifelse(test = (file_type == "coverage"),
                             yes = str_c("track type=bigWig name=\"", bw_name, "\" bigDataUrl=\"", URL, "\""),
                             no = str_c("track type=bam name=\"", bam_name, "\" bigDataUrl=\"", URL, "\""))) %>%
  dplyr::select(-c(bw_name, bam_name)) %>%
  tidyr::unite(scale_type, scaled, file_type, sep = ".") 

# create table with all possible columns
tracks_placeholder <- 
  tibble(sample_id = rep(unique(tracks_tbl$sample_id), each = 3),
         scale_type = rep(c("raw.coverage", "RPM_scaled.coverage", "raw.individual_reads"), length(unique(tracks_tbl$sample_id))))

# join placeholder with table
tracks_tbl_tidy <- 
  tracks_tbl %>% 
  dplyr::select(-experiment) %>% 
  dplyr::right_join(., tracks_placeholder, by = c("sample_id", "scale_type")) %>% 
  tidyr::spread(key = scale_type, value = URL) %>% 
  dplyr::mutate(experiment = experiment)

# combine and save
stats_and_tracks <-
  right_join(stats_tbl, tracks_tbl_tidy, by = "sample_id") %>%
  dplyr::select(experiment, sample_id, raw.coverage, RPM_scaled.coverage, raw.individual_reads, everything()) %T>%
  readr::write_csv(., path = file.path(outpath, str_c("log.", experiment, ".stats_and_tracks.csv")))


### plot
# data for plot
plot_tb <- 
  stats_and_tracks %>% 
  dplyr::select(sample_id, raw_input:other) %>% 
  dplyr::mutate(sample_id = sample_id %>% 
                  str_remove_all(., "^s_|\\.SE$|\\.PE$") %>% 
                  str_replace_all(., "_", " "))

# total mapped reads
mapped_total <- 
  plot_tb %>% 
  dplyr::select(sample_id, raw_input:unmapped_total) %>% 
  dplyr::mutate(mapped_total = round((mapped_total / raw_input), 2) * 100, 
                unmapped_total = round((unmapped_total / raw_input), 2) * 100) %>% 
  dplyr::select(-raw_input) %>% 
  tidyr::pivot_longer(cols = -sample_id, names_to = "value", values_to = "count") %>% 
  dplyr::mutate(value = factor(value, levels = c("unmapped_total", "mapped_total")))

# plot stacked barplot
stacked_bar <- 
  ggplot() +
  geom_bar(data = mapped_total, aes(x = sample_id, y = count, fill = value), position = "stack", stat = "identity") +
  scale_fill_viridis_d() + 
  guides(fill = F) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) +
  theme(legend.title = element_blank())

# save plot
ggsave(plot = stacked_bar, 
       filename = file.path(outpath, "test.png"), width = 10, height = 10)



