library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(readr)
library(stringr)

# library size data.frame
logs_filenames <- file.path(c("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_1cell.WE/s_1cell.WELog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE/s_2cell.WELog.final.out", 
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_2cell.WE_DNAm/s_2cell.WE_DNAmLog.final.out",
                              "/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_4cell.WE/s_4cell.WELog.final.out"))

library_size_df <- 
  data.frame(sapply(X = 1:length(logs_filenames), function(X) as.integer(read.delim(logs_filenames[X], header = F, stringsAsFactors = F)[8, 2]))) %>% 
  set_colnames("library_size") %>% 
  mutate(sample_name = c("s_GV", "s_1C", "s_2C", "s_2Ca", "s_4C"), 
         library_size = library_size / 10^6) %>% 
  dplyr::select(2:1)
  
# set tile width
tile_width <- 10000

# normalize to fpkm
all_coverage_fpkm <- 
  all_coverage_counts_df_merged %>% 
  set_colnames(c("pos", library_size_df$sample_name)) %>% 
  tidyr::gather(key = sample_name, value = fpkm, -pos) %>% 
  left_join(library_size_df, by = "sample_name") %>% 
  mutate(fpkm = fpkm / (library_size * (tile_width / 1000)), 
         sample_name = factor(sample_name, levels = library_size_df$sample_name)) %>% 
  dplyr::select(sample_name, pos, fpkm) %>% 
  tidyr::spread(key = sample_name, value = fpkm)

# divide to upstream/downstream
# add bin column 
all_coverage_fpkm_binned <- 
  rbind(all_coverage_fpkm %>% #upstream
          filter(pos >= -150000 & pos < 0) %>%
          mutate(bin = rev(gl(ceiling(n() / tile_width), tile_width, n())), 
                 stream = "upstream", 
                 pos = abs(pos)),
        # all_coverage_fpkm %>% #element
        #   filter(pos >= 0 & pos <= element_width) %>%
        #   mutate(bin = 1),  
        all_coverage_fpkm %>% #downstream
          filter(pos > element_width & pos <= (150000 + element_width)) %>%
          mutate(bin = gl(ceiling(n() / tile_width), tile_width, n()), 
                 stream = "downstream", 
                 pos = pos - element_width)) %>% 
  mutate(bin = as.integer(bin))

# take individual stage data and do 
stage_name <- "s_2C"

wilcox_results <- 
  all_coverage_fpkm_binned %>% 
  dplyr::select_("pos", "bin", "stream", "fpkm" = stage_name) %>% 
  reshape2::dcast(pos + bin ~ stream, value.var = "fpkm") %>% 
  group_by(bin) %>% 
  # bootstrap(m = 1000, by_group = T) %>%
  do(tidy(wilcox.test(.$upstream, .$downstream, alternative = "less", paired = T))) %>% 
  ungroup(wilcox_results) %>% 
  as.data.frame()
