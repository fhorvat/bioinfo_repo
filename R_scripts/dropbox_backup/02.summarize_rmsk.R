### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes/rmsk_counts/Analysis")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(xlsx)
library(ggthemes)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set base path
base_path <- "/common/WORK/fhorvat/Projekti/Svoboda/maternal_transcriptomes"

######################################################## READ DATA

######################################################## MAIN CODE
## set experiments
experiment_list <- c("chinese_hamster.CHOK1GS_HDv1", "cow.bosTau8", "golden_hamster.mesAur1", "mouse_B6.mm10", 
                     "mouse_cast.CAST_EiJ", "mouse_PWD.PWK_Phj", "pig.susScr11", "rabbit.oryCun2", "rat.rn6")

# set animal names
experiment_table <- 
  tibble(experiment_id = experiment_list) %>% 
  dplyr::mutate(animal = str_remove(experiment_id, "\\..*") %>% str_replace_all(., "_", " "), 
                animal = str_c(toupper(str_sub(animal, 1, 1)), str_sub(animal, 2, nchar(animal)), sep = ""), 
                animal = factor(animal, levels = c("Mouse B6", "Mouse cast", "Mouse PWD", "Rat", 
                                                   "Chinese hamster", "Golden hamster", "Rabbit", "Cow", "Pig")))

### loop through experiments, clean rmsk RDS, join all experiments to one list
rmsk_list <- purrr::map(experiment_list, function(experiment){
  
  ### set paths
  # repeatMasker RDS path
  rmsk_rds_path <- file.path(base_path, "rmsk_counts", experiment)
  
  # mapped path
  mapped_path <- file.path(base_path, "Mapped", experiment)
  
  # classified reads path
  read_class_path <- list.files(path = rmsk_rds_path, pattern = "rmsk.read_class.*RDS", full.names = T)
  
  # stats&tracks path
  stats_path <- list.files(mapped_path, pattern = "*stats_and_tracks.csv", full.names = T)
  
  
  ### read data
  # read class reads table, gets different mismatch counts
  read_class_list <- purrr::map(read_class_path, readRDS)
  
  # read stats and tracks table 
  library_size_df <- 
    readr::read_csv(file = stats_path) %>% 
    dplyr::select(sample_id, mapped_total)
  
  
  ### clean repeatMasker .RDS
  read_class_clean <- 
    purrr::map(c("repClass", "repFamily", "repName"), function(repeat_annotation){
      
      purrr::map(read_class_list, function(sample) sample[[repeat_annotation]]) %>% 
        dplyr::bind_rows(.) %>% 
        dplyr::mutate(sample_id = str_remove(sample_id, "\\.total")) %>% 
        dplyr::mutate_at(vars(starts_with("rep")), funs(str_remove(., "\\?"))) %>% 
        dplyr::group_by_("sample_id", repeat_annotation) %>% 
        dplyr::summarise(count = sum(count)) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::left_join(., library_size_df, by = "sample_id") %>% 
        dplyr::mutate(total_fraction = (count / mapped_total)) %>% 
        dplyr::group_by(sample_id) %>% 
        dplyr::mutate(repeat_fraction = (count / sum(count))) %>% 
        dplyr::group_by_(repeat_annotation) %>% 
        dplyr::summarise(repeat_fraction = mean(repeat_fraction), 
                         total_fraction = mean(total_fraction)) %>% 
        dplyr::arrange(desc(repeat_fraction)) %>% 
        dplyr::mutate(experiment_id = experiment)
      
    }) %>% 
    magrittr::set_names(c("repClass", "repFamily", "repName"))
  
})

# join all repeat annotations from different experiments to one data.frame
rmsk_annotation <- purrr::map(c("repClass", "repFamily", "repName"), function(repeat_annotation){
  
  purrr::map(rmsk_list, function(x) x[[repeat_annotation]]) %>% 
    dplyr::bind_rows(.)
  
}) %>% 
  magrittr::set_names(c("repClass", "repFamily", "repName"))

# repFamily
repfamily_plot_df <- 
  rmsk_annotation$repFamily %>% 
  group_by(grp) %>%
  top_n(n = 5, wt = x)

  dplyr::mutate(repClass = replace(repClass, repClass == "Simple_repeat", "simple repeat"), 
                repeat_fraction = repeat_fraction * 100, 
                repClass = ifelse(repClass %in% repClass_categories, repClass, "other repeat")) %>% 
  dplyr::group_by(repClass, experiment_id) %>% 
  dplyr::summarise(repeat_fraction = sum(repeat_fraction)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., experiment_table, by = "experiment_id") %>% 
  dplyr::mutate(repClass = factor(repClass, levels = c(repClass_categories, "other repeat")), 
                repeat_fraction = round(repeat_fraction, 1))

### visualize
# choose categories to visualize, set everything else to "other repeat"
repClass_categories <- c("SINE", "LINE", "LTR", "simple repeat", "rRNA")

# repClass
repclass_plot_df <- 
  rmsk_annotation$repClass %>% 
  dplyr::mutate(repClass = replace(repClass, repClass == "Simple_repeat", "simple repeat"), 
                repeat_fraction = repeat_fraction * 100, 
                repClass = ifelse(repClass %in% repClass_categories, repClass, "other repeat")) %>% 
  dplyr::group_by(repClass, experiment_id) %>% 
  dplyr::summarise(repeat_fraction = sum(repeat_fraction)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::left_join(., experiment_table, by = "experiment_id") %>% 
  dplyr::mutate(repClass = factor(repClass, levels = c(repClass_categories, "other repeat")), 
                repeat_fraction = round(repeat_fraction, 1))

# plot
ggplot(data = repclass_plot_df, mapping = aes(repClass, repeat_fraction, fill = animal)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 1)) +
  geom_text(mapping = aes(label = str_c(repeat_fraction, "%"), y = repeat_fraction, color = animal), 
            angle = 90, size = 5, hjust = 1, position = position_dodge(width = 0.9)) +
  scale_color_manual(values = c(rep("white", 7), rep("black", 2))) +
  scale_fill_viridis_d() + 
  guides(fill = guide_legend(nrow = 1), 
         color = F) +
  theme_tufte() +
  theme(axis.text.x = element_text(angle = 0, size = 15, vjust = 1),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(), 
        axis.title = element_blank(),
        legend.text = element_text(size = 15), 
        legend.position = "bottom") +
  ggsave(file.path(outpath, "repClass_percent.facet.png"), width = 15, height = 10)

