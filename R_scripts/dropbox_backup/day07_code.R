### INFO: Advent of code 2017, day01 (http://adventofcode.com/2017/day/1)
### DATE: 18. 12. 2017 
### AUTHOR: Filip Horvat

################################################################################### SCRIPT PARAMS
rm(list = ls()); gc()

################################################################################### WORKING DIRECTORY
setwd("C:/Users/Filip/Dropbox/Bioinfo/other_projects/test/adventOfCode/2017")

################################################################################### LIBRARIES
# data shaping
library(magrittr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)

################################################################################### SOURCE FILES

################################################################################### FUNCTIONS

################################################################################### PATH VARIABLES
# set in and out path
inpath <- getwd()

################################################################################### TABLES

################################################################################### MAIN CODE
### part 1 
input <- readr::read_lines(file = "day07_input.txt")

bottom <- 
  tibble(program = input[str_detect(input, "->")]) %>% 
  tidyr::separate(col = program, into = c("program", "weight"), sep = "\\(") %>% 
  dplyr::mutate(weight = str_replace(weight, "\\)", "")) %>% 
  tidyr::separate(col = weight, into = c("weight", "subprogram"), sep = "->")

subprograms <- str_split(bottom$subprogram, pattern = ",") %>% unlist() %>% stringr::str_trim()
programs <- stringr::str_trim(bottom$program)
root_name <- programs[!(programs %in% subprograms)]

### part 2
# get programs and their weigths
programs_weights <-
  tibble(subprogram = input) %>%
  dplyr::mutate(subprogram = str_replace(subprogram, "->.*", "") %>% str_trim(.), 
                program_name = str_replace(subprogram, "\\(.*\\)", "") %>% str_trim(.)) 

# get programs and their relations
programs_relations <- 
  tibble(program = input) %>% 
  tidyr::separate(col = program, into = c("program", "subprogram"), sep = "->") %>% 
  dplyr::mutate(program = str_trim(program), 
                subprogram = str_trim(subprogram))

# get root
tree <- 
  programs_relations %>% 
  dplyr::filter(str_detect(program, root_name)) %>%   
  tidyr::separate(col = subprogram, into = str_c("sub", 1:(str_count(.$subprogram, ",") + 1)), sep = ", ") %>% 
  tidyr::gather(value, subprogram, -program) %>% 
  dplyr::select(-value) %>% 
  dplyr::left_join(., programs_weights, by = c("subprogram" = "program_name")) %>% 
  dplyr::select(program, subprogram = subprogram.y) %>% 
  data.table::setnames(., old = ncol(.), "last")
  
# fill tree in loop
while(TRUE){
  
  # create temp tree
  tree_temp <- 
    tree %>% 
    dplyr::select(program = last) %>% 
    dplyr::filter(!is.na(program)) %>% 
    left_join(., programs_relations, by = "program") %>%
    dplyr::filter(!is.na(subprogram))

  # break if temp_tree is empty
  if(nrow(tree_temp) == 0){
    break
  }
  
  tree_temp %<>% 
    tidyr::separate(col = subprogram, into = str_c("sub", 1:(max(str_count(.$subprogram, ",")) + 1)), sep = ", ") %>% 
    tidyr::gather(value, subprogram, -program) %>% 
    dplyr::select(-value) %>% 
    dplyr::left_join(., programs_weights, by = c("subprogram" = "program_name")) %>% 
    dplyr::select(program, subprogram = subprogram.y)
  
  # join with original tree
  tree <- 
    left_join(tree, tree_temp, by = c("last" = "program")) %>% 
    data.table::setnames(., old = 2:ncol(.), c(str_c("sub", 1:(ncol(.) - 2)), "last"))
  
}

tree <-  
  unique(tree) %>%
  data.table::setnames(., str_c("sub", ncol(.):1)) %>% 
  dplyr::mutate_all(.funs = funs(replace(., is.na(.), "bla (0)")))

### 
while(TRUE){
  
  tree_temp <-
    tree %>% 
    dplyr::select(sub2, sub1) %>%
    dplyr::mutate(weight1 = str_extract(sub1, "\\(.*\\)") %>% str_replace_all(., "\\(|\\)", "") %>% as.integer(.),
                  weight2 = str_extract(sub2, "\\(.*\\)") %>% str_replace_all(., "\\(|\\)", "") %>% as.integer(.)) %>% 
    dplyr::group_by(sub2) %>% 
    dplyr::summarise(weight1 = sum(weight1),
                     weight2 = unique(weight2)) %>% 
    dplyr::mutate(weight = weight1 + weight2) %>% 
    dplyr::select(sub2, weight) %>% 
    dplyr::filter(sub2 != "bla (0)")
  
  tree_short <-
    dplyr::left_join(tree, tree_temp, by = "sub2") %>%
    dplyr::filter(weight != 0) %>%
    dplyr::group_by(sub3) %>%
    dplyr::summarise(unique_weight = str_c(unique(weight), collapse = "|")) %>%
    dplyr::filter(str_count(unique_weight, "\\|") > 0)
  
  if(nrow(tree_short) > 0){
    break
  }
  
  tree <- 
    dplyr::left_join(tree, tree_temp, by = "sub2") %>% 
    dplyr::mutate(weight = str_c(str_replace(sub2, " \\(.*\\)", ""), " (", weight, ")")) %>% 
    dplyr::select(-sub2, -sub1) %>% 
    data.table::setnames(., str_c("sub", ncol(.):1)) %>% 
    unique(.) %>% 
    dplyr::mutate_all(.funs = funs(replace(., is.na(.), "bla (0)")))
}

tree_problem <-
  tree %>%
  dplyr::filter(sub3 == tree_short$sub3) %>%
  dplyr::select(sub2, sub1) %>%
  dplyr::mutate(weight1 = str_extract(sub1, "\\(.*\\)") %>% str_replace_all(., "\\(|\\)", "") %>% as.integer(.),
                weight2 = str_extract(sub2, "\\(.*\\)") %>% str_replace_all(., "\\(|\\)", "") %>% as.integer(.)) %>%
  dplyr::group_by(sub2) %>%
  dplyr::summarise(weight1 = sum(weight1),
                   weight2 = unique(weight2)) %>%
  dplyr::mutate(weight_total = weight1 + weight2) %>%
  dplyr::select(sub2, weight2, weight_total) %>%
  dplyr::filter(sub2 != "bla (0)")

weight_wrong <- 
  table(tree_problem$weight_total) %>% 
  as.tibble(.) %>% 
  dplyr::arrange(desc(n)) %>% 
  slice(2) %$%
  Var1 %>% 
  as.numeric(.)

weight_difference <- 
  table(tree_problem$weight_total) %>% 
  as.tibble(.) %>% 
  dplyr::arrange(desc(n)) %>% 
  slice(1) %$%
  Var1 %>% 
  as.numeric(.) %>% 
  magrittr::subtract(weight_wrong, .)

### solution
tree_problem %>% 
  dplyr::filter(weight_total == weight_wrong) %$%
  weight2 %>%
  magrittr::subtract(., weight_difference)

