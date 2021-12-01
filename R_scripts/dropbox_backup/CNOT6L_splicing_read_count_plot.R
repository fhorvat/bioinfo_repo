library(dplyr)
library(readr)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tibble)

library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(BiocParallel)
library(doMC)
library(genomation)

### replicating figure 6D from EMBO-J paper - comparison of unspliced/spliced read pair ratios per cell stage
# UNSPLICED 
# - read counts where one end maps to intron/exon junction or entirely in intron and the
#   other end maps to the adjacent exon were labeled as ‘unspliced’
# SPLICED 
# - one end maps either to the splice site and covers two adjacent exons or 
#   with each end mapping to separate, adjacent exons

# some code from vfranke:
# /home/members/vfranke/Projects/PSvoboda_Oocyte_Transcriptome/Scripts/GenomeActivation/SplicingAnalysis/SplicingAnalysis_Scripts

################################################################## functions
### loads an RData file and assigns it to the name variable
Assigner = function(`.path`, `.name`){
  
  if(!is.character(`.path`) | !is.character(`.name`)){
    stop('Both arguments should be characters!')
  }
  
  load(`.path`)
  assign(`.name`, get(ls()[1]), parent.frame())
  
}

################################################################## parameters and paths 
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads")

#  register the multicore parallel backend with the foreach package 
registerDoMC(21)

# set paths
inpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/R_objects"
outpath <- "/common/WORK/fhorvat/Projekti/Svoboda/Analyses/CNOT6L/spliced_reads/results"

################################################################## read data
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
  mutate(ID = str_replace_all(sample_path, "^/.*/|.bam", "")) %>%
  rowwise() %>%
  mutate(library_size = as.integer(read_delim(logs_path, delim = "\t", col_names = F)[8, 2])) %>% 
  dplyr::select(ID, library_size, sample_path) 

################################################################## 
rfiles <- list.files(inpath, full.names = TRUE, pattern = 'RData')

# load RData files
l <- foreach(i = 1:length(rfiles)) %dopar% {
  
  print(i)
  Assigner(rfiles[i], "a")
  
  return(a)
}

# set names and order
names(l) <- str_replace(basename(rfiles), ".SpliceDonAcc.RData", "")
# l <- l[ord]

# extract second data.frame from each element of the list
sl <- lapply(l, function(x) x[[2]])
sl <- lapply(names(sl), function(x){
  a <- sl[[x]]
  a$samp <- x
  return(a)
  })

# bind all data.frames in list to one data.frame
sl <- rbindlist(sl)
sl$samp <- factor(sl$samp, levels = names(l))
# setnames(sl, 3:8, str_replace(colnames(sl)[3:8], "genes.", ""))
sl$ex[sl$int] <- FALSE

# get library size
tot.counts <- 
  as.integer(library_size_df$library_size) %>% 
  set_names(., library_size_df$ID)
tot.counts <- tot.counts[names(tot.counts) %in% names(l)]
tot.norm <- round(1e7 / tot.counts, 2)

################################################################## plot
d <- sl[, list(ex = sum(ex), int = sum(int)), by = samp]
m <- melt(d)
m$norm <- m$value * tot.norm[m$samp]
m.rat <- m[, norm[variable == 'int'] / norm[variable == 'ex'], by = samp]

ggplot(m.rat, aes(x = samp, y = V1, fill = "int.ex.ratio")) +
  geom_bar(stat = "identity", color = "black", fill = "black") +
  theme_bw() + 
  # theme(panel.border = element_rect(size = .1, colour = "black"),
  #       axis.ticks = element_line(size = 1),
  #       axis.text = element_text(size = 20),
  #       axis.text.x = element_text(angle = 90),
  #       axis.title = element_text(size = 20),
  #       panel.grid = element_blank(),
  #       legend.key.size = unit(1, "cm"),
  #       legend.text = element_text(size = 15),
  #       legend.title = element_text(size = 15),
  #       legend.position = "none") +
  ylab("Unspliced / Spliced") +
  xlab("") +
  ggsave(filename = file.path(outpath, "Splicing.ReadCount_test1.pdf"))

