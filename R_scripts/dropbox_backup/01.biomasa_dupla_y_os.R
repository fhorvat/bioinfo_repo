### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/tmp/test/Antonija")

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
# read files
data <- read.csv(file="Biomasa.csv")
data$Number_of_cell_per_L <- data$Number_of_cell_per_L / 10000000

p <- ggplot(data, aes(x = sample, group=1))
p <- p + geom_line(aes(y = Biomass_mgL.1, colour = "Biomass_mgL/L"))

p <- p + geom_line(aes(y = Number_of_cell_per_L, colour = "Number_of_cells/L")) 
p <- p + scale_y_continuous(sec.axis = sec_axis(~ .* 10000000, name = "Number_of_cells/L", labels = scales::comma))


p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Biomass_mg/L",
              x = "Sample",
              colour = "Parameter")+
  scale_x_discrete(limits=c("S06","S07","S08","S09", "S10", "S11", "S12", "S01", "S02", "S03"))
p <- p + theme(legend.position = c(0.8, 0.9))
p  <- p+ theme_classic()


ggsave("p.png", p, width = 8, height =6 )
