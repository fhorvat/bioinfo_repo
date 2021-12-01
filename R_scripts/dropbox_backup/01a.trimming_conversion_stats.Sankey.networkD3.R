### INFO: 
### DATE: Fri Oct 09 09:41:14 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_bisulfite.test_run/Data/Raw/Cleaned.bbduk_adapter")

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

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# list trimming logs paths
trim_logs_path <- list.files(inpath, ".*\\.trim\\.log", full.names = T)

# conversion read number path
conversion_log <- file.path(inpath, "converted_reads", "fastq_files.counts.txt")

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
                read_in_pair = str_extract(sample_id, "1$|2$") %>% str_c("read_", .), 
                sample_id = str_remove(sample_id, "_[1,2]$")) %>% 
  tidyr::pivot_wider(., id_cols = sample_id, names_from = c("tmp", "read_in_pair"), values_from = "count") %>% 
  dplyr::select(-trimmed_read_2) %>% 
  dplyr::select(sample_id, from_trimmed = trimmed_read_1, converted_read_1, converted_read_2)

######################################################## MAIN CODE
# join tables 
bisulfite_stats <- 
  trim_logs_tb %>% 
  dplyr::left_join(., conversion_log_tb, by = "sample_id") %>% 
  dplyr::mutate(trimmed_reads = from_trimmed * 2, 
                converted_reads = converted_read_1 + converted_read_2, 
                removed_conversion_filtering = trimmed_reads - converted_reads) %>% 
  dplyr::select(sample_id, 
                input_reads = Input, 
                removed_adapter_trimming = Removed, pass_adapter_trimming = trimmed_reads, 
                removed_conversion_filtering, pass_conversion_filtering = converted_reads)

### plot as Sankey plot
# create template links table
links_template <- tibble(source = c("input_reads", "input_reads", 
                                    "pass_adapter_trimming", "pass_adapter_trimming"),
                         target = c("pass_adapter_trimming", "removed_adapter_trimming", 
                                    "pass_conversion_filtering", "removed_conversion_filtering"))

# choose sample
sample <- "s_hamster_GV_Mov10l1_HET_So787_r1.PE"

# create table
links_tb <- 
  bisulfite_stats %>% 
  dplyr::filter(sample_id == sample) %>% 
  dplyr::mutate_at(.vars = vars(!matches("sample_id")), .funs = as.numeric) %>% 
  tidyr::pivot_longer(., cols = -sample_id, names_to = "target", values_to = "count") %>% 
  dplyr::select(-sample_id) %>% 
  dplyr::left_join(links_template, ., by = "target")

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name = c(as.character(links_tb$source), as.character(links_tb$target)) %>% unique(.))

# add groups
nodes$group <- as.factor(c("default"))

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links_tb$IDsource <- match(links_tb$source, nodes$name) - 1 
links_tb$IDtarget <- match(links_tb$target, nodes$name) - 1
links_tb$group <- as.factor(c("a", "b", "a", "b"))

# prepare color scale
my_color <- 'd3.scaleOrdinal() .domain(["a", "b", "default"]) .range(["steelblue", "gray", "gray"])'

# Make the Network
p <- sankeyNetwork(Links = links_tb, 
                   Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "count", NodeID = "name", 
                   colourScale = my_color, 
                   NodeGroup = "group", LinkGroup = "group",
                   sinksRight = FALSE, 
                   fontSize = 50, 
                   nodeWidth = 20,
                   nodePadding = 30, 
                   height = 600, 
                   width = 2000)

# save as network
p %>% 
  htmlwidgets::prependContent(htmltools::tags$h1("Title")) %>% 
  saveNetwork(file = 'title_Mis.html')

# save the widget
saveWidget(p, file = file.path(outpath, "sankeyBasic1.html"))

