### INFO: 
### DATE: Thu Apr 18 17:45:13 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/transcriptome_assemblies/cow.bosTau8/trinity.extended_merge.bbnorm_normalized")

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
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# assembly stats file
stats_path <- list.files(inpath, "*Trinity.stats", full.names = T)

# assembly fasta index file
index_path <- list.files(inpath, "*Trinity.fasta.fai", full.names = T)

# busco short summary path - mammalia
busco_mammalia_path <- list.files(inpath, "short_summary.*.Trinity.busco.mammalia.txt", full.names = T, recursive = T)

# busco short summary path - metazoa
busco_metazoa_path <- list.files(inpath, "short_summary.*.Trinity.busco.metazoa.txt", full.names = T, recursive = T)

# minimap2 alignment - reads to transcriptome
read_counts_path <- list.files(file.path(inpath, "minimap2_align.reads_to_transcriptome"), "*read_counts.txt", full.names = T)

# minimap2 alignment - transcripts to genome
read_counts_path <- list.files(file.path(inpath, "minimap2_align.transcriptome_to_genome"), "*read_counts.txt", full.names = T)

######################################################## READ DATA
# read assembly stats file
stats_tb <- readr::read_delim(stats_path, 
                              delim = ":", 
                              skip = 5, 
                              col_name = c("category", "count"))

# read assembly fasta index file
index_tb <- readr::read_delim(index_path, 
                              delim = "\t", 
                              col_names = c("name", "length", "offset", "linebases", "linewidth"))

# read busco mammalia report
busco_mammalia_tb <- readr::read_delim(busco_mammalia_path, 
                                       delim = "\t", 
                                       skip = 7,
                                       col_names = c("tmp", "count", "category"))

# read busco metazoa report
busco_metazoa_tb <- readr::read_delim(busco_metazoa_path, 
                                      delim = "\t", 
                                      skip = 7,
                                      col_names = c("tmp", "count", "category"))

# read counts - reads to transcriptome
counts_reads_transcriptome_tb <- purrr::map(read_counts_path, function(path){
  
  readr::read_delim(path, 
                    delim = "\t", 
                    col_names = "count") %>% 
    mutate(pairing = str_extract(path, "SE|PE"))
  
}) %>% 
  bind_rows(.)

# read count - transcripts to genome
counts_transcripts_genome_tb <- readr::read_delim(read_counts_path, 
                                                  delim = "\t", 
                                                  col_names = "count")

######################################################## MAIN CODE
# get name of the sample
sample_name <- 
  basename(stats_path) %>% 
  str_remove(., "\\.Trinity\\.stats")

### clean tables
# basic stats
stats_tb %<>% 
  filter(., !is.na(count)) %>% 
  mutate(count = count %>% str_remove(., "\t") %>% str_trim(.) %>% as.integer(.), 
         category = category %>% str_remove(., "\t") %>% str_replace_all(., " ", "_") %>% str_remove_all(., "Percent_|Contig_|'") %>% tolower(.), 
         description = c(rep("general", 3), rep("all_transcripts", 8), rep("longest_isoform", 8)), 
         category = description %>% 
           str_c(., ".", category) %>% 
           str_remove(., "^general\\.") %>% 
           str_replace(., "gc", "GC_perc") %>% 
           str_replace(., "\\.n", ".N"))

## Busco
# mammalia
busco_mammalia_tb %<>% 
  select(-tmp) %>% 
  filter(!is.na(category)) %>% 
  mutate(perc = ((as.integer(count) / 4104) * 100), 
         species = str_c(sample_name, ".mammalia")) %>% 
  filter(category != "Total BUSCO groups searched")

# metazoa
busco_metazoa_tb %<>% 
  select(-tmp) %>% 
  filter(!is.na(category)) %>% 
  mutate(perc = ((as.integer(count) / 978) * 100), 
         species = str_c(sample_name, ".metazoa")) %>% 
  filter(category != "Total BUSCO groups searched")

# clean counts of reads mapping to transcriptome
counts_reads_transcriptome_tb_clean <- 
  counts_reads_transcriptome_tb %>% 
  mutate(mapping = rep(c("mapped", "total"), 2), 
         category = str_c(pairing, mapping, sep = ".")) %>% 
  select(-c(pairing, mapping)) %>% 
  tidyr::spread(key = category, value = count) %>% 
  mutate(mapped = PE.mapped + SE.mapped, 
         total = PE.total + SE.total, 
         perc = round((mapped / total), 3) * 100)


# # plot BUSCO
# ggplot() +
#   geom_bar(data = busco_mammalia_tb %>% filter(category != "Complete BUSCOs (C)"), aes(y = perc, x = species, fill = category), stat = "identity", width = 0.75) +
#   # coord_flip() +
#   theme_gray(base_size = 8) +
#   scale_y_continuous(labels = c("0","20","40","60","80","100"), 
#                      breaks = c(0, 20, 40, 60, 80, 100)) +
#   scale_fill_manual(values = c("#56B4E9", "#3492C7", "#F0E442", "#F04442"),
#                     labels = c("Complete (C) and single-copy (S)",
#                                "Complete (C) and duplicated (D)",
#                                "Fragmented (F)",
#                                "Missing (M)")) +
#   ggtitle("BUSCO Assessment Results - mammalia") +
#   xlab("") +
#   ylab("\n%BUSCOs") +
#   theme(plot.title = element_text(family = "sans", colour = "black", face = "bold"), 
#         # legend.position = "top", 
#         legend.title = element_blank(), 
#         legend.text = element_text(family = "sans", size = rel(1.2) * 1), 
#         panel.background = element_rect(color = "#FFFFFF", fill = "white"), 
#         panel.grid.minor = element_blank(), 
#         panel.grid.major = element_blank(), 
#         axis.text.y = element_text(family = "sans", colour = "black", size = rel(1.66) * 1),
#         axis.text.x = element_text(family = "sans", colour = "black", size = rel(1.66) * 1), 
#         axis.line = element_line(size = 1*1, colour = "black"), 
#         axis.ticks.y = element_line(colour = "white", size = 0), 
#         axis.ticks.x = element_line(colour = "#222222"), 
#         axis.ticks.length = unit(0.4, "cm"), 
#         axis.title.x = element_text(family = "sans", size = rel(1.2) * 1)) +
#   guides(fill = guide_legend(override.aes = list(colour = NULL))) +
#   ggsave(file.path(outpath, "test.png"), width = 15, height = 10)

