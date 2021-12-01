### INFO:
### DATE: Tue Jan 29 20:37:35 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
setwd(".")

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

library(GenomicRanges)
library(GenomicAlignments)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# plot
collapseReads <- function(bam){
  
  ### count overlaps
  # all
  counts_all <- countOverlaps(line1_tiled, bam, ignore.strand = T)
  
  # sense
  counts_sense <- countOverlaps(line1_tiled, bam, ignore.strand = F)
  
  # antisense
  counts_antisense <- counts_all - counts_sense
  
  ### sum across positions in LINE1
  # sum sense
  counts_sense %<>%
    matrix(., nrow = 50) %>%
    rowSums(.)
  
  # sum antisense
  counts_antisense %<>%
    matrix(., nrow = 50) %>%
    rowSums(.)
  
  ### create table and plot
  # table
  plot_df <- tibble(pos = c(1:50, 50:1),
                    counts = c(counts_sense, - counts_antisense),
                    direction = rep(c("sense", "antisense"), each = 50),
                    bam_name = bam_name,
                    experiment_name = experiment_name)
  
  # return
  return(plot_df)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# list of LINE1 full length elements path
line1_coords_path <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/Documentation/LINE1_full_length.20180517.ZJM.tidy.csv"

# bam path
bam_path <- "%BAM_PATH"

# get names
bam_name <- basename(bam_path) %>% str_remove(., "\\.[SE|PE].*")
experiment_name <- dirname(bam_path) %>% basename(.)

######################################################## READ DATA
# read Zoe's list of LINE1 full length elements
line1_coords <- read_csv(line1_coords_path)

######################################################## MAIN CODE
### prepare data
# form GRanges
line1_gr <- GRanges(line1_coords)

# tile each LINE1 to 50 tiles (= 2% of length)
line1_tiled <- tile(line1_gr, n = 50)
names(line1_tiled) <- line1_gr$id
line1_tiled <- unlist(line1_tiled)

### small RNA data
# read bam, get each read only once
bam_gr <-
  readGAlignmentsList(file = bam_path,
                      use.names = TRUE,
                      param = ScanBamParam(which = reduce(line1_gr))) %>%
  unlist(.) %>%
  sortSeqlevels(.) %>%
  sort(.) %>%
  .[!duplicated(names(.))]

# 21-23 (siRNA) reads, weighted
collapseReads(bam = bam_gr[width(bam_gr) >= 21 & width(bam_gr) <= 23]) %>%
  readr::write_csv(., file.path(outpath, "results", str_c("full_LINE1.cummulative", experiment_name, bam_name, "21to23nt.weighted", "csv", sep = ".")))

# 24-31 (siRNA) reads, weighted
collapseReads(bam = bam_gr[width(bam_gr) >= 24 & width(bam_gr) <= 31]) %>%
  readr::write_csv(., file.path(outpath, "results", str_c("full_LINE1.cummulative", experiment_name, bam_name, "24to31nt.weighted", "csv", sep = ".")))

