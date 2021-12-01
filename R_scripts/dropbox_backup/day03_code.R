### INFO: Advent of Code 2018, day 01
### DATE: Sat Dec 01 13:44:36 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()

######################################################## WORKING DIRECTORY
setwd("C:/Users/fhorvat/Dropbox/Bioinfo/other_projects/test/adventOfCode/2018")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

######################################################## READ DATA
# read data
input <- readr::read_lines(file = file.path(inpath, "day03_input.txt"))

######################################################## MAIN CODE
#### PART 1 ####
# add to tibble
input_tb <- 
  tibble(raw_input = input) %>% 
  tidyr::separate(raw_input, into = c("id", "temp", "coords", "dim"), sep = " ", remove = F) %>% 
  tidyr::separate(dim, into = c("x_dim", "y_dim"), sep = "x") %>% 
  tidyr::separate(coords, into = c("x_start", "y_end"), sep = ",") %>% 
  dplyr::select(-temp) %>% 
  dplyr::mutate(id = str_remove(id, "^#"),
                y_end = str_remove(y_end, ":$")) %>% 
  dplyr::mutate_at(vars(matches("x|y")), funs(as.integer(.))) %>% 
  dplyr::mutate(y_end = 1000 - y_end) %>% 
  dplyr::mutate(x_end = x_start + x_dim, 
                y_start = y_end - y_dim) %>% 
  dplyr::select(id, x_start, x_end, y_start, y_end)

# plot
ggplot() +
  geom_rect(data = input_tb, mapping = aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
            fill = "red", alpha = 0.5, color = "black") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 100)) +
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 100))

### fill matrix
# list input tb
input_list <- 
  input_tb %>% 
  dplyr::select(-id) %>% 
  dplyr::mutate(x_start = x_start + 1, 
                y_start = 1000 - y_start, 
                y_end = 1000 - y_end + 1) %>% 
  t() %>% 
  as.data.frame(.) %>% 
  as.list(.)

# matrix
mat <- matrix(0, nrow = 1000, ncol = 1000)

# fill matrix
for(n in 1:length(input_list)){
  
  mat[input_list[[n]][3]:input_list[[n]][4], input_list[[n]][1]:input_list[[n]][2]] <- mat[input_list[[n]][3]:input_list[[n]][4], input_list[[n]][1]:input_list[[n]][2]] + 1
  
}

sum(mat > 1)

 
#### PART 2 ####
input <- readr::read_lines(file = file.path(inpath, "day03_input.txt"))

# clean input
input_tb_clean <- 
  tibble(raw_input = input) %>% 
  tidyr::separate(raw_input, into = c("id", "temp", "coords", "dim"), sep = " ", remove = F) %>% 
  tidyr::separate(dim, into = c("x_dim", "y_dim"), sep = "x") %>% 
  tidyr::separate(coords, into = c("x_start", "y_end"), sep = ",") %>% 
  dplyr::select(-temp) %>% 
  dplyr::mutate(id = str_remove(id, "^#"),
                y_end = str_remove(y_end, ":$")) %>% 
  dplyr::mutate_at(vars(matches("x|y")), funs(as.integer(.)))

# calculate start and end coordinates for each ID
input_tb <- 
  input_tb_clean %>% 
  dplyr::mutate(x_start = x_start + 1, 
                x_end = x_start + x_dim - 1, 
                y_start = y_end + y_dim, 
                y_end = y_end + 1, 
                total_area = x_dim * y_dim) %>% 
  dplyr::select(id, x_start, x_end, y_start, y_end, x_dim, y_dim, total_area)

# list input tb
input_list <- 
  input_tb %>% 
  dplyr::select(x_start, x_end, y_start, y_end) %>% 
  t() %>% 
  as.data.frame(.) %>% 
  as.list(.)

# matrix
mat <- matrix("x", nrow = 1000, ncol = 1000)

# fill matrix
for(n in 1:length(input_list)){
  
  mat[input_list[[n]][3]:input_list[[n]][4], input_list[[n]][1]:input_list[[n]][2]] <- 
    str_c(mat[input_list[[n]][3]:input_list[[n]][4], input_list[[n]][1]:input_list[[n]][2]], n) 
  
}

# remove x, count, compare
mat_tb <- 
  table(mat) %>% 
  as.tibble(.) %>% 
  set_colnames(., c("id", "freq")) %>% 
  dplyr::mutate(id = str_remove(id, "x")) %>% 
  dplyr::right_join(., input_tb %>% dplyr::select(id, total_area), by = "id") %>% 
  dplyr::filter(freq == total_area) %$%
  id
