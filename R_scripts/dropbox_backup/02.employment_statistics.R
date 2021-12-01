### INFO: NGSchool 2019 survey
### DATE: Wed Feb 20 08:05:40 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/test/NGSchool")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)
library(openxlsx)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# table path
table_path <- file.path(inpath, "national_M2017_dl.xlsx")

######################################################## READ DATA
# read table
employment_tb <- 
  read.xlsx(table_path) %>% 
  as.tibble(.)

######################################################## MAIN CODE
# C2 - There is a weak negative correlation (cor = -0.12, p-value = 9.1e-04) 
# between the number of employees and their median salary. I.e. the more people employed in the profession, 
# the less they seem to earn. You can identify both, the most exclusive, profitable and popular, 
# low-paying occupation by calculating the ratio of the salary to the number of employees. 
# Provide the difference in annual median salary between the two opposite extremes of such ratio. 
# Restrict your analysis to the detailed groups.The answer is an integer; type it in without any white spaces, commas etc.

# tidy and filter table
extreme_salaries <- 
  employment_tb %>% 
  select(occ_title = OCC_TITLE, occ_group = OCC_GROUP, total_employees = TOT_EMP, median_annual_salary = A_MEDIAN) %>% 
  filter(occ_group == "detailed", 
         !median_annual_salary %in% c("*", "#")) %>% 
  mutate(median_annual_salary = as.numeric(median_annual_salary), 
         ratio_salary_employees = median_annual_salary / total_employees) %>% 
  filter((ratio_salary_employees == max(.$ratio_salary_employees)) | (ratio_salary_employees == min(.$ratio_salary_employees))) %>% 
  arrange(median_annual_salary) %$%
  median_annual_salary

# answer
extreme_salaries[2] - extreme_salaries[1]


# C4 - Let's assume you are about to enter the job market and you want to make informed decision about your career. 
# At this time, you are solely interested in the financial aspect of your future job. You are pretty confident you will 
# achieve a lot in any field but also want to be sure that even if things won't go according to plan you end up earning 
# better than in another job. You want to choose a job which has the highest 90th percentile of hourly wage but which 
# H_PCT25 is in the upper 10% of all 25th percentile hourly wages. To make it simple you are constraining your analysis 
# to the major groups - you'll look into more details some other time. Please, provide the OCC_TITLE of the occupation 
# you come up during your analysis.

# tidy and filter
best_job <- 
  employment_tb %>% 
  select(occ_title = OCC_TITLE, occ_group = OCC_GROUP, h_pct90 = H_PCT90, h_pct25 = H_PCT25) %>% 
  filter(occ_group == "major") %>% 
  mutate_at(vars(starts_with("h_")), funs(as.integer(.))) %>% 
  filter_at(vars(starts_with("h_")), all_vars(!is.na(.))) %>% 
  arrange(desc(h_pct25)) %>% 
  slice(1:(round(nrow(.) * 0.1))) %>% 
  filter(h_pct90 == max(h_pct90)) %$% 
  occ_title

