### INFO: 
### DATE: Mon Oct 28 18:43:00 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
# 
######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/PhD/algorithms_and_programming/2019_10_28/homework")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# breeding table path
breeding_path <- file.path(inpath, "Breeding_Performance_Lnc1_190815.csv")

######################################################## READ DATA
# read genes info
breeding_tb <- fread(breeding_path)

######################################################## MAIN CODE
# a)
breeding_tb <- fread(breeding_path)

breeding_tb_tidy <- 
  breeding_tb %>% 
  as_tibble(.) %>% 
  dplyr::rename_all(str_replace_all, " ", "_") %>%
# or any of the following:  
# dplyr::rename_at(.vars = vars(contains("Age")), .funs = ~ str_replace_all(., " ", "_"))
# dplyr::rename_all(.funs = ~ str_replace_all(., " ", "_"))
# magrittr::set_colnames(., str_replace_all(colnames(.), " ", "_"))
# dplyr::rename(`Age of mother in days` = Age_of_mother_in_days, 
#               `Age of mother In weeks` = Age_of_mother_In_weeks) 
  dplyr::rename(father_genotype = male, mother_genotype = female) %>% 
  dplyr::select(-c(Timeline:D.O.L))

# b) also ok if they put ratio to new column with mutate and then filter 
breeding_tb_tidy %>% 
  dplyr::filter(Total > Dead, 
                Female >= (Male * 2.5)) %>% 
  dplyr::arrange(desc(Age_of_mother_in_days))

# b) 
breeding_tb_tidy %>% 
  group_by(Cage) %>% 
  dplyr::summarise(KO = sum(`-/-`)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::filter(KO == max(KO, na.rm = T)) %$% 
  Cage
   
# c) if they didn't notice that there is one row where Total = 0, still give them points. If only 1-2 of them noticed, maybe give extra points?
breeding_tb_tidy %>% 
  dplyr::filter(Dead == 0, Total > 0) %>% 
  dplyr::mutate(mating_age = Age_of_mother_in_days - 21) %$% 
  mating_age %>% 
  mean(.)

# d)
breeding_tb_tidy %>% 
  dplyr::group_by(Cage) %>% 
  dplyr::summarise(male = sum(Male), 
                   female = sum(Female))

# e)
breeding_tb_tidy %>% 
  dplyr::filter(mother_genotype == father_genotype) %>% 
  dplyr::group_by(mother_genotype, father_genotype) %>% 
  dplyr::summarise(dead = sum(Dead))

# f) 
breeding_tb_tidy %>% 
  dplyr::mutate(age_category = ifelse(Age_of_mother_In_weeks < 36, "young", "old")) %>% 
  dplyr::group_by(age_category) %>% 
  dplyr::summarise_at(vars("Male", "Female"), funs(mean)) %>% 
  pivot_longer(data = ., 
               cols = -age_category, 
               names_to = "pup_sex",
               values_to = "mean")
