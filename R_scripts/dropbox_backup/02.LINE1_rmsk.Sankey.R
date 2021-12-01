### INFO: 
### DATE: Thu Oct 03 16:18:30 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/other/lab_meetings/Prague/2019_10_04.LINE1")

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

library(networkD3)
library(htmlwidgets)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# whole joined repeatMasker path
rmsk_path <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2/rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz"

# get annotated LINE1 path
line1_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Dicer_Mili_KO/Analysis/RNAi_piRNA_paper/filter_LINE1/LINE1.4000nt_plus.annotated.csv"

######################################################## READ DATA
# read RepeatMasker
rmsk_tb <- 
  readr::read_delim(rmsk_path, delim = "\t") %>% 
  dplyr::mutate(strand = "*") %>% 
  GRanges(.) %>% 
  as_tibble(.)

# read LINE1's
line1_tb <- readr::read_csv(line1_path)

######################################################## MAIN CODE
# create links
links <- 
  tibble(source = c(
    # "RepeatMasker all insertions", "RepeatMasker all insertions", 
    "LINE1 insertions", "LINE1 insertions", 
    "whole", "whole", 
    "full insertion", "full insertion", 
    "no exon overlap", "no exon overlap"), 
    target = c(
      # "LINE1 insertions", "other insertions", 
      "whole", "interupted", 
      "full insertion", "shorter or longer", 
      "no exon overlap", "overlaps exon", 
      "has both ORFs", "ORFs too short"),
    value = c(
      # 636811, 2442097,
      503158, 133653, 
      15295, 487863, 
      14975, 320, 
      1578, 13397)) %>% 
  as.data.frame(.)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name = c(as.character(links$source), as.character(links$target)) %>% unique(.))

# add groups
nodes$group <- as.factor(c("default"))

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name) - 1 
links$IDtarget <- match(links$target, nodes$name) - 1
links$group <- as.factor(c("a", "b", "a", "b", "a", "b", "a", "b"))

# prepare color scale
my_color <- 'd3.scaleOrdinal() .domain(["a", "b", "default"]) .range(["steelblue", "gray", "gray"])'

# Make the Network
p <- sankeyNetwork(Links = links, 
                   Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   colourScale = my_color, 
                   NodeGroup = "group", LinkGroup = "group",
                   sinksRight = FALSE, 
                   fontSize = 20, 
                   nodeWidth = 20)

# save the widget
saveWidget(p, file = file.path(outpath, "sankeyBasic1.html"))

