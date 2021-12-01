#!/common/WORK/vfranke/bin/R/R-3.1.0/bin/Rscript
### INFO: R Script
### DATE: 05.08.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads")

# {1} LIBRARIES
lib.path='/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/Fugaku/vfranke_documentation'
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)
library(genomation)
library(ggplot2)
library(readr)

#/{1} LIBRARIES


# {2} CODE
# {{1}} FUNCTIONS

#/{{1}} FUNCTIONS


# {{2}} INPUT VARIABLES 

# {{{1}}} PATH VARIABLES
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/R_objects"
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/results"

#/{{{1}}} PATH VARIABLES

# {{{2}}} SCRIPT PARAMS
registerDoMC(21)


# {{3}} MAIN

# CNOT6L library size
library_size_df <- 
  tibble(sample_path = list.files("/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                  pattern = "*.bam$",
                                  recursive = T, 
                                  full.names = T), 
         logs_path = list.files(path = "/common/WORK/fhorvat/Projekti/Svoboda/CNOT6L/Data/Mapped/STAR_mm10", 
                                pattern = "*Log.final.out", 
                                recursive = T, 
                                full.names = T)) %>%
  dplyr::mutate(ID = str_replace_all(sample_path, "^/.*/|_.*", "")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size, sample_path)

# get library size
tot.counts <- 
  as.integer(library_size_df$library_size) %>% 
  set_names(., str_replace(string = library_size_df$ID, pattern = "_.*", replacement = ""))
tot.counts <- tot.counts[names(tot.counts) %in% names(l)]
tot.norm <- round(1e7 / tot.counts, 2)


# read files
rfiles = list.files(inpath, full.names=TRUE, pattern='RData')
l = foreach(i = 1:length(rfiles))%dopar%{
  
  print(i)
  Assigner(rfiles[i], 'a')
  return(a)
}
names(l) = str_replace(basename(rfiles), "_.*", "")
# l = l[ord]

sl = lapply(l, function(x)x[[2]])
sl = lapply(names(sl), function(x){a=sl[[x]];a$samp=x;a})
sl = rbindlist(sl)
sl$ex[sl$int] = FALSE

# plot
d = sl[,list(ex=sum(ex), int=sum(int)),by=samp]
m = melt(d[!str_detect(d$samp,'PA'),])
m$norm = m$value*tot.norm[m$samp]
m.rat = m[,norm[variable=='int']/norm[variable=='ex'], by=samp]
setnames(m.rat,2,'ratio')
ggplot(m.rat, aes(x=samp, y=ratio)) +
  geom_bar(stat="identity") +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10)) + 
  ggtitle('MII removed, Intron/Exon ratio') +
  ggsave("test_plot_2.pdf", width = 9.93, height = 5.86)


