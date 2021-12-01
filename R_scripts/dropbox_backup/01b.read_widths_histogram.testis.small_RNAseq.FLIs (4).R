### INFO: 
### DATE: Thu Jul 30 18:45:22 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/IAP/library_composition")

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
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/expression/hamster_testis_Mov10l.small_RNAseq"

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# mapped path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_testis_Mov10l.small_RNAseq"
documentation_path <- file.path(base_path, "Data/Documentation")
sample_table_path <- list.files(documentation_path, "\\.sampleTable\\.csv", full.names = T)
mapped_path <- file.path(base_path, "Data/Mapped/STAR_Siomi.multimappers")
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

# gene info path
gene_info_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP"
gene_info_path <- file.path(gene_info_path, "IAP.FLI_elements.csv")

# expression files path
expression_path <- file.path(inpath, "expression_sum.sense.RDS_files")
expression_path <- list.files(expression_path, ".*\\.expression\\.RDS$", full.names = T)

######################################################## READ DATA
# read expression values
expression_tb_full <- purrr::map(expression_path, function(path){
  
  # get sample name
  sample_name <- path %>% basename(.) %>% str_remove(., "\\.expression\\.RDS")
  
  # read table
  tmp_tb <- 
    path %>% 
    readRDS(.) %>% 
    dplyr::mutate(sample_id = sample_name)
  
  # return
  return(tmp_tb)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(repName = "IAP_full")

# read sample table
sample_tb <- readr::read_csv(sample_table_path) 

# read library sizes
library_size_tb <- readr::read_delim(library_size_path, delim = "\t", col_names = c("sample_id", "library_size")) 

# filter and summarize library size
library_size_tb %<>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_") %>% 
  dplyr::group_by(genotype_age) %>% 
  dplyr::summarise(library_size = sum(library_size))

# filter sample table
sample_tb %<>% 
  dplyr::select(sample_id, genotype, age) %>% 
  tidyr::unite(col = "genotype_age", genotype, age, sep = "_") 

######################################################## MAIN CODE
### full length guys
# set name
plot_name <- "05.intact_guys"

# filter table
expression_tb <- expression_tb_full


### get expression 
# add category to expression table
count_tb <- 
  expression_tb %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt")) %>% 
  dplyr::group_by(sample_id, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(!str_detect(sample_id, "HET")) %>% 
  dplyr::left_join(., sample_tb, by = "sample_id") %>% 
  dplyr::group_by(genotype_age, read_width) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  # tidyr::pivot_wider(id_cols = c("genotype_age"), names_from = c(read_width), values_from = count) %>% 
  dplyr::left_join(., library_size_tb, by = c("genotype_age")) %>% 
  dplyr::mutate(library_size = library_size / 1e6, 
                RPM = count / library_size) %>% 
  dplyr::mutate(genotype_age = factor(genotype_age, levels = c("Mov10l_WT_13dpp", "Mov10l_KO_13dpp", "Mov10l_WT_21dpp", "Mov10l_KO_21dpp")))

# save as table
count_tb %>% 
  tidyr::pivot_wider(., id_cols = "genotype_age", names_from = read_width, values_from = RPM, names_prefix = "r.") %>% 
  dplyr::arrange(genotype_age) %T>%
  readr::write_csv(str_c(plot_name, "testis", "small_RNAseq", "RPM", "read_length.histogram.csv", sep = "."))

# plot
read_length_plot <- 
  ggplot() +
  geom_histogram(data = count_tb, 
                 mapping = aes(x = read_width, y = RPM, fill = genotype_age), width = 0.9, stat = "identity") +
  facet_grid(genotype_age ~ ., scales = "free") +
  scale_x_discrete(limits = 19:32, breaks = 19:32, labels = as.character(19:32)) +
  guides(fill = guide_legend(reverse = F)) +
  ggtitle(plot_name) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  theme(strip.text.y = element_text(size = 6, colour = "black"),
        strip.placement = "outside")

# save
ggsave(plot = read_length_plot, filename = file.path(outpath, 
                                                     str_c(plot_name, "testis", "small_RNAseq", "RPM", "read_length.histogram.png", sep = ".")), 
       width = 10, height = 7.5)
