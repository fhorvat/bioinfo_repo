library(readr)
library(dplyr)
library(magrittr)
library(stringr)

# paths
input_dir <- "/common/WORK/fhorvat/Projekti/Svoboda/lncRNA_KO_2016/Analysis/variant_calling/GATK/01_STAR_1st_map"
sjout_path <- file.path(input_dir, "all_samples_SJ.out.tab")

# filter 
sjout <- 
  read_delim(file = sjout_path, delim = "\t", col_names = F) %>% 
  set_colnames(c("chr", "start", "end", "strand", "intron_motif", "splice_annotation", "unique_reads", "multimapping_reads", "max_splice_overhang")) %>% 
  filter(splice_annotation == 0, 
         chr != "chrM") %T>%
  readr::write_delim(., path = file.path(input_dir, "all_samples_filtered_SJ.out.tab"), col_names = F, delim = "\t")
