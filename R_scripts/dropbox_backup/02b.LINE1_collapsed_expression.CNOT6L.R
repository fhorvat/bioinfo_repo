### INFO: 
### DATE: Tue Jan 29 20:37:35 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression")

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
plotCummulativeReads <- function(bam, bam_path, title){
  
  # get names
  bam_name <- basename(bam_path) %>% str_remove(., "\\.[SE|PE].*")
  experiment_name <- str_remove(bam_path, "\\/Data.*") %>% basename(.) %>% str_remove(., "_.*")
  
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
                    direction = rep(c("sense", "antisense"), each = 50))
  
  # plot
  ggplot() +
    geom_rect(data = plot_df, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = counts, fill = direction)) +
    geom_hline(yintercept = 0, color = "black") +
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggsave(filename = file.path(outpath, "results", str_c("full_LINE1.cummulative", experiment_name, bam_name, title, "png", sep = ".")),
           width = 15,
           height = 10)
  
}

######################################################## PATH VARIABLES
# set main inpath
inpath <- getwd()

# set main outpath
outpath <- getwd()

# list of LINE1 full length elements path
line1_coords_path <- file.path(inpath, "Documentation", "LINE1_full_length.20180517.ZJM.tidy.csv")

# bam path
bam_path <- "/common/WORK/fhorvat/Projekti/Svoboda/accessory_data_sets/small_RNAseq/Yang_2016_SciAdv_GSE83581/Data/Mapped/STAR_mm10.filtered/s_oocyte_r1.SE.19to32nt.bam" 

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
                      param = ScanBamParam(which = reduce(line1_gr), 
                                           flag = scanBamFlag(isMinusStrand = NA))) %>% 
  unlist(.) %>% 
  sortSeqlevels(.) %>%
  sort(.) %>%
  .[!duplicated(names(.))]

# 21-23 (siRNA) reads, weighted
plotCummulativeReads(bam = bam_gr[width(bam_gr) >= 21 & width(bam_gr) <= 23], 
                     bam_path = bam_path, 
                     title = "21to23nt.weighted")

# 24-31 (piRNA reads), weighted 
plotCummulativeReads(bam = bam_gr[width(bam_gr) >= 24 & width(bam_gr) <= 31], 
                     bam_path = bam_path, 
                     title = "24to31nt.weighted")


### CNOT6L
# path
bam_path_CNOT6L <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/retrotransposon_expression/datasets/CNOT6L/Mapped/1_full_LINE1_reads/s_GV_WT_r1.PE.total.full_LINE1.bam"

# read CNOT6L GV WT bam
bam_gr_CNOT6L <-
  readGAlignmentsList(file = bam_path_CNOT6L, use.names = TRUE)%>%
  unlist(.) %>% 
  sortSeqlevels(.) %>% 
  sort(.) %>% 
  grglist(.) %>%
  .[!duplicated(names(.))] %>%
  unlist(.)

# plot all length CNOT6L reads
plotCummulativeReads(bam = bam_gr_CNOT6L, bam_path = bam_path_CNOT6L, title = "CNOT6L.weighted")
