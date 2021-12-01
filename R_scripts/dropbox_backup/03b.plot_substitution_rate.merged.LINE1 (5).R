### INFO: 
### DATE: Wed Jul 15 10:34:22 2020
### AUTHOR: fhorvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LINE1/substitution_rate/random_sampled_LINE1s")

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

# table with LINE1s info path
line1_tb_path <- file.path(inpath, "LINE1s.random_200_per_repName.200805.csv")

# substitution rate .RDS files path
sub_rate_path <- file.path(inpath, "sub_rate.RDS_files")
sub_rate_path <- list.files(sub_rate_path, pattern = "^sub_rate\\..*\\.RDS$", full.names = T)

######################################################## READ DATA
# read table with LINE1s info
line1_tb <- readr::read_csv(line1_tb_path)

# read substitution rate files, bind to one table
sub_rate_tb_all <- purrr::map(sub_rate_path, function(path){
  
  # get name
  fasta_name <- path %>% basename(.) %>% str_remove_all(., "^sub_rate\\.|\\.RDS$")
  
  # read from .RDS
  sub_rate_raw <- readRDS(path)
  
  # add to tibble
  sub_rate_raw_tb <- tibble(sub_rate = sub_rate_raw, 
                            repName = fasta_name)
  
}) %>% 
  bind_rows(.)

######################################################## MAIN CODE
### merge some repName
# create classes list
comboClasses <- list("HAL1" = c("HAL1", "HAL1b", "HAL1M8", "HAL1ME", "MusHAL1"),
                     "L1Lx" = c("L1Lx_I", "L1Lx_II", "L1Lx_III", "L1Lx_IV"), 
                     "L1M"  = c("L1M", "L1M1", "L1M2", "L1M2a", "L1M2b", "L1M2c", "L1M3", "L1M3a", "L1M3b", 
                                "L1M3c", "L1M3d", "L1M3de", "L1M3e", "L1M3f", "L1M4", "L1M4a1", "L1M4a2", 
                                "L1M4b", "L1M4c", "L1M5", "L1M6", "L1M6B", "L1M7", "L1M8"),
                     "L1MA" = c("L1MA10", "L1MA4", "L1MA4A", "L1MA5", "L1MA5A", "L1MA6", "L1MA7", "L1MA8", "L1MA9"),
                     "L1MB" = c("L1MB1", "L1MB2", "L1MB3", "L1MB4", "L1MB5", "L1MB7", "L1MB8"),
                     "L1MC" = c("L1MC", "L1MC1", "L1MC2", "L1MC3", "L1MC4", "L1MC4a", "L1MC5", "L1MC5a", "L1MCa", "L1MCb", "L1MCc"),
                     "L1MD" = c("L1MD", "L1MD1", "L1MD2", "L1MD3", "L1MDa", "L1MDb"),
                     "L1MdA" = c("L1MdA_I", "L1MdA_II", "L1MdA_III", "L1MdA_IV", "L1MdA_V", "L1MdA_VI", "L1MdA_VII"),
                     "L1MdF" = c("L1MdFanc_I", "L1MdFanc_II", "L1MdF_I", "L1MdF_II", "L1MdF_III", "L1MdF_IV", "L1MdF_V"), 
                     "L1MdGf" = c("L1MdGf_I", "L1MdGf_II"), 
                     "L1MdMus" = c("L1MdMus_I", "L1MdMus_II"), 
                     "L1MdN" = c("L1MdN_I"),
                     "L1MdTf" = c("L1MdTf_I", "L1MdTf_II", "L1MdTf_III"), 
                     "L1MdV" = c("L1MdV_I", "L1MdV_II", "L1MdV_III"), 
                     "L1ME" = c("L1MEa", "L1MEb", "L1MEc", "L1MEd", "L1MEf", "L1MEg", "L1MEg1", "L1MEg2", "L1MEh", "L1MEi", "L1MEj"), 
                     "L1ME1" = c("L1ME1"), 
                     "L1ME2" = c("L1ME2", "L1ME2z"),
                     "L1ME3" = c("L1ME3", "L1ME3A", "L1ME3B", "L1ME3C", "L1ME3Cz", "L1ME3D", "L1ME3E", "L1ME3F", "L1ME3G"), 
                     "L1ME4" = c("L1ME4a", "L1ME4b", "L1ME4c"), 
                     "L1ME5" = c("L1ME5"), 
                     "L1_Rod" = c("L1_Rod"), 
                     "L1_Mur" = c("L1_Mur1", "L1_Mur2", "L1_Mur3"), 
                     "L1_Mus" = c("L1_Mus3"),
                     "L1P" = c("L1P5", "L1PB4"), 
                     "L1VL1" = c("L1VL1"),
                     "Lx" = c("Lx"), 
                     "Lx2" = c("Lx2", "Lx2A", "Lx2A1", "Lx2B", "Lx2B2"), 
                     "Lx3" = c("Lx3A", "Lx3B", "Lx3C", "Lx3_Mus"), 
                     "Lx4" = c("Lx4A", "Lx4B"), 
                     "Lx5" = c("Lx5", "Lx5b", "Lx5c"), 
                     "Lx6" = c("Lx6"), 
                     "Lx7" = c("Lx7"),
                     "Lx8" = c("Lx8", "Lx8b"), 
                     "Lx9" = "Lx9", 
                     "Lx10" = "Lx10")

# list to tibble
classes_tb <- purrr::map(names(comboClasses), function(ltr_name){
  
  # get one LTR type in tibble
  comboClasses[[ltr_name]] %>% 
    tibble(repName = ., repSubfamily = ltr_name)
  
}) %>%
  dplyr::bind_rows(.)

# save
readr::write_csv(classes_tb, file.path(outpath, str_c("LINE1s.random_200_per_repSubfamily.200805", "classes.csv", sep = ".")))

### list of retrotransposons to plot
# prepare table for plot 
sub_rate_tb <- 
  sub_rate_tb_all %>% 
  dplyr::right_join(., classes_tb, by = "repName") %>%
  dplyr::group_by(repSubfamily) %>% 
  dplyr::slice_sample(n = 200) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(repSubfamily = factor(repSubfamily, levels = unique(classes_tb$repSubfamily))) %>%
  dplyr::arrange(repSubfamily)

# data statistics
sub_rate_stats <- 
  sub_rate_tb %>% 
  group_by(repSubfamily) %>% 
  summarize(median = median(sub_rate), 
            sd = sd(sub_rate), 
            length = length(sub_rate)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

# labels for plot
labels_tb <- 
  sub_rate_tb %>% 
  dplyr::group_by(repSubfamily) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count > 0)

# plot
sub_rate_boxplot <- 
  ggplot() +
  stat_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), geom = "errorbar", size = 1.5) +
  geom_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), outlier.colour = NULL, outlier.shape = NA, size = 1.5) +
  geom_jitter(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), alpha = 0.4, 
              colour = "black", size = 1.5, width = 0.1, height = 0, show.legend = F) +
  scale_x_discrete(labels = str_c(labels_tb$repSubfamily, 
                                  " (",
                                  labels_tb$count, 
                                  ")"), 
                   drop = TRUE) + 
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())

# save plot
ggsave(plot = sub_rate_boxplot, 
       filename = file.path(outpath, str_c("LINE1s.random_200_per_repSubfamily.200805", "boxplot.png", sep = ".")), 
       width = 30, 
       height = 10)

