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

######################################################## READ DATA

######################################################## MAIN CODE
### merge some repName
# create classes list from golden hamster
comboClasses.golden_hamster <- list("HAL1" = c("HAL1", "HAL1b", "HAL1M8", "HAL1ME", "MusHAL1"),
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
                                    "L1Md_T" = c("L1MdTf_I", "L1MdTf_II", "L1MdTf_III"), 
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
                                    "L1VL" = c("L1VL1"),
                                    "Lx" = c("Lx"), 
                                    "Lx2" = c("Lx2", "Lx2A", "Lx2A1", "Lx2B", "Lx2B2"), 
                                    "Lx3" = c("Lx3A", "Lx3B", "Lx3C", "Lx3_Mus"), 
                                    "Lx4" = c("Lx4A", "Lx4B"), 
                                    "Lx5" = c("Lx5", "Lx5b", "Lx5c"), 
                                    "Lx6" = c("Lx6"), 
                                    "Lx7" = c("Lx7"),
                                    "Lx8" = c("Lx8", "Lx8b"), 
                                    "Lx9" = c("Lx9"), 
                                    "Lx10" = c("Lx10"))

# LINE1s in mouse but but not annotated in golden hamster
comboClasses.mouse <- list("HAL1" = c("HAL1_SS"),
                           "L1M" = c("L1M2a1"), 
                           "L1MdA" = c("L1Md_A"), 
                           "L1MdF" = c("L1Md_F", "L1Md_F2",  "L1Md_F3"), 
                           "L1MdGf" = c("L1Md_Gf"),
                           "L1Md_T" = c("L1Md_T"), 
                           "L1_Mus" = c("L1_Mm", "L1_Mus1", "L1_Mus2", "L1_Mus4"), 
                           "L1VL" = c("L1VL2", "L1VL4"))

# LINE1s in rat but not but not annotated in golden hamster
comboClasses.rat <- list("HAL1" = c("RNHAL1"),
                         "L1M" = c("L1M2a1"), 
                         "L1VL" = c("L1VL4a", "L1VL4"), 
                         "Lx3" = c("Lx3"))

# LINE1s in chinese hamster but not annotated in golden hamster
comboClasses.chinese_hamster <- list("L1M" = c("L1M2a1"), 
                                     "L1_Rod" = c("L1_Rod2"))


### lists to one tibble
classes_tb_all <- purrr::map(list(comboClasses.golden_hamster, comboClasses.mouse, comboClasses.rat, comboClasses.chinese_hamster), function(comboClasses){
  
  # join all repNames in one list
  classes_tb <- purrr::map(names(comboClasses), function(ltr_name){
    
    # get one LTR type in tibble
    class_ltr <- 
      comboClasses[[ltr_name]] %>% 
      tibble(repName = ., repSubfamily = ltr_name)
    
    # return
    return(class_ltr)
    
  }) %>%
    dplyr::bind_rows(.)
  
  # return
  return(classes_tb)
  
}) %>% 
  dplyr::bind_rows() %>% 
  unique(.)


# save table
readr::write_csv(classes_tb_all, file.path(outpath, "LINE1s.random_200_per_repSubfamily.classes.200807.csv"))

