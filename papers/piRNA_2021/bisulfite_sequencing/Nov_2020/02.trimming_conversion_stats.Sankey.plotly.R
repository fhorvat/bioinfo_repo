### INFO: 
### DATE: Fri Oct 09 09:41:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Raw/Cleaned")

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

library(networkD3)
library(htmlwidgets)
library(plotly)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get ggplot colors
gg_color_hue <- function(n){
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list trimming logs paths
trim_logs_path <- list.files(inpath, ".*\\.trim\\.log", full.names = T)

# conversion read number path
conversion_log <- file.path(inpath, "converted_reads", "fastq_files.counts.txt")

# mapping stats log
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped"
mapping_log <- list.files(mapped_path, "\\.stats_and_tracks\\.csv", full.names = T, recursive = T)
mapping_log <- mapping_log[str_detect(mapping_log, "converted_unstranded_PE")]

######################################################## READ DATA
# read trimming logs
trim_logs_tb <- purrr::map(trim_logs_path, function(path){
  
  # read lines
  trim_log <- 
    readr::read_lines(path) %>%
    .[str_detect(., "in1=|Input:|KTrimmed:|Total Removed:")] %>% 
    tibble(raw = .) %>% 
    dplyr::mutate(clean = raw %>% str_remove_all(., ".*in1=|\t| reads.*|:| in2\\=.*|Total |,") %>% str_squish(.) %>% basename(.)) %>% 
    dplyr::select(clean) %>% 
    unique(.) %>% 
    dplyr::mutate(clean = ifelse(!str_detect(clean, " "), str_c("sample_id ", clean), clean)) %>% 
    tidyr::separate(clean, into = c("category", "value"), sep = " ") %>% 
    tidyr::pivot_wider(., names_from = "category", values_from = "value")
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "(?<=\\.[P,S]E).*"))

# read conversion log
conversion_log_tb <- 
  readr::read_delim(conversion_log, delim = "\t", col_names = c("sample_id", "count")) %>% 
  dplyr::mutate(tmp = str_extract(sample_id, "converted$"), 
                tmp = replace(tmp, is.na(tmp), "trimmed"), 
                sample_id = str_remove(sample_id, "\\.converted$"), 
                read_in_pair = str_extract(sample_id, "1$|2$|s$") %>% str_c("read_", .), 
                sample_id = str_remove(sample_id, "_[1,2,s]$")) %>% 
  tidyr::pivot_wider(., id_cols = sample_id, names_from = c("tmp", "read_in_pair"), values_from = "count")
  # dplyr::select(-trimmed_read_2) %>% 
  # dplyr::select(sample_id, from_trimmed = trimmed_read_1, converted_read_1, converted_read_2, converted_read_s)

# # mapping logs
# mapping_log_tb <- purrr::map(mapping_log, function(path){
#   
#   # read report log
#   report_tb <- 
#     readr::read_csv(path) %>%
#     dplyr::mutate(mapped_reads = `Aligned Reads` + `Ambiguously Aligned Reads`, 
#                   unmapped_reads = `Unaligned Reads`) %>% 
#     dplyr::select(sample_id, mapped_reads, unmapped_reads)
#   
# }) %>% 
#   bind_rows(.) %>% 
#   dplyr::filter(!is.na(mapped_reads)) %>% 
#   dplyr::mutate(sample_id = str_remove(sample_id, "(?<=\\.[P,S]E).*")) %>% 
#   dplyr::group_by(sample_id) %>% 
#   dplyr::summarise(mapped_reads = sum(mapped_reads),
#                    unmapped_reads = sum(unmapped_reads)) %>% 
#   dplyr::ungroup(.)

######################################################## MAIN CODE
# join tables 
bisulfite_stats <- 
  trim_logs_tb %>% 
  dplyr::left_join(., conversion_log_tb, by = "sample_id") %>% 
  # dplyr::left_join(., mapping_log_tb, by = "sample_id") %>% 
  dplyr::mutate(trimmed_reads = from_trimmed * 2, 
                converted_reads = converted_read_1 + converted_read_2 + converted_read_s, 
                removed_conversion_filtering = trimmed_reads - converted_reads) %>% 
  dplyr::select(sample_id, 
                input_reads = Input, 
                removed_adapter_trimming = Removed, pass_adapter_trimming = trimmed_reads, 
                removed_conversion_filtering, pass_conversion_filtering = converted_reads, 
                pass_mapping = mapped_reads, removed_during_mapping = unmapped_reads)

### plot as Sankey plot
# create template links table
links_template <- tibble(source = c("input_reads", "input_reads", 
                                    "pass_adapter_trimming", "pass_adapter_trimming", 
                                    "pass_conversion_filtering", "pass_conversion_filtering"),
                         target = c("pass_adapter_trimming", "removed_adapter_trimming", 
                                    "pass_conversion_filtering", "removed_conversion_filtering", 
                                    "pass_mapping", "removed_during_mapping"))

# loop through samples
purrr::map(unique(bisulfite_stats$sample_id), function(sample_name){
  
  # create table
  links_tb <- 
    bisulfite_stats %>% 
    dplyr::filter(sample_id == sample_name) %>% 
    dplyr::mutate_at(.vars = vars(!matches("sample_id")), .funs = as.numeric) %>% 
    dplyr::rename_at(.vars = vars(!matches("sample_id")), ~str_c(., ".count")) %>% 
    dplyr::mutate(input_reads.percentage = 100, 
                  removed_adapter_trimming.percentage = round(100 * (removed_adapter_trimming.count / input_reads.count), 2), 
                  pass_adapter_trimming.percentage = round(100 * (pass_adapter_trimming.count / input_reads.count), 2), 
                  removed_conversion_filtering.percentage = round(100 * (removed_conversion_filtering.count / input_reads.count), 2), 
                  pass_conversion_filtering.percentage = round(100 * (pass_conversion_filtering.count / input_reads.count), 2), 
                  removed_during_mapping.percentage = round(100 * (removed_during_mapping.count / input_reads.count), 2), 
                  pass_mapping.percentage = round(100 * (pass_mapping.count / input_reads.count), 2)) %>% 
    tidyr::pivot_longer(., cols = -sample_id, names_to = c("target", "value_type"), names_sep = "\\.", values_to = "value") %>% 
    tidyr::pivot_wider(., id_cols = c(sample_id, target), names_from = "value_type", values_from = "value") %>% 
    dplyr::select(-sample_id) %>% 
    dplyr::left_join(links_template, ., by = "target")
  
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- 
    tibble(name = c(as.character(links_tb$source), as.character(links_tb$target)) %>% unique(.)) %>% 
    dplyr::mutate(node_color = ifelse(str_detect(name, "pass|input"), "steelblue", "gray")) %>% 
    dplyr::left_join(., links_tb %>% dplyr::select(target, percentage), by = c("name" = "target")) %>% 
    dplyr::mutate(percentage = replace(percentage, is.na(percentage), 100)) %>% 
    dplyr::mutate(x_pos = c(0.1, 0.35, 0.6, 0.35, 0.6, 0.9, 0.9), 
                  y_pos = 0.5, 
                  y_pos = ifelse(str_detect(name, "pass"), (0.5 * (percentage / 100)), (1 - (0.5 * (percentage / 100)))))
  
  # connection must be provided using id, not using real name like in the links dataframe
  links_tb$IDsource <- match(links_tb$source, nodes$name) - 1 
  links_tb$IDtarget <- match(links_tb$target, nodes$name) - 1
  
  # prepare color scale
  plot_colors <- 
    gg_color_hue(length(nodes$name)) %>% 
    col2rgb(., alpha = F) %>% 
    as.vector() %>% 
    split(., ceiling(seq_along(.) / 3)) %>% 
    purrr::map(., function(x) str_c(x, collapse = ",") %>% str_c("rgba(", ., ",0.8)")) %>% 
    unlist(.) %>% 
    unname(.)
  
  # Make the Network
  fig <- plot_ly(type = "sankey",
                 arrangement = "perpendicular",
                 domain = list(x =  c(0, 1), 
                               y =  c(0, 1)),
                 orientation = "h",
                 valueformat = ".0f",
                 valuesuffix = "\nreads",
                 
                 node = list(label = str_replace_all(nodes$name, "_", " ") %>% str_c(., " (", nodes$percentage, "%)"),
                             color = nodes$node_color,
                             x = nodes$x_pos, 
                             y = nodes$y_pos, 
                             pad = 5,
                             thickness = 15,
                             line = list(color = "black",
                                         width = 0.5), 
                             hovertemplate = '%{value:,}<extra></extra>'),
                 
                 link = list(source = links_tb$IDsource,
                             target = links_tb$IDtarget,
                             value =  links_tb$count, 
                             hovertemplate = '%{value:,}<extra></extra>'))
  
  # add title and other layout parameters
  fig <- 
    fig %>% 
    layout(title = sample_name,
           font = list(size = 10),
           xaxis = list(showgrid = F, zeroline = F),
           yaxis = list(showgrid = F, zeroline = F))
  
  # save as html
  saveWidget(fig, file = file.path(outpath, str_c(sample_name, "read_stats", "Sankey_diagram.html", sep = ".")))
  
  # return
  return(sample_name)
  
})


