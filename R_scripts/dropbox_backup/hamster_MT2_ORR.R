library("dplyr")
library("readr")
library("ggplot2")
library("GenomicRanges")

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/hamster")

# reading hamster repeatMasker table
mesAur1_rmsk <- 
  read_delim("/common/DB/vfranke/Base/GenomeAnnotation/mesAur1/Refseq/mesAur1.RepeatMasker.Refseq.txt",
             delim = "\t", 
             col_names = c("seqnames", "start", "end", "strand", "element_name", "element_class")) %>%
  mutate(strand = replace(strand, strand == "C", "*"))

# filtering LTRs
mesAur1_LTRs <- 
  mesAur1_rmsk %>%
  filter(grepl("LTR", element_class)) %>%
  mutate(width = end - start + 1) %>%
  mutate(strand = replace(strand, strand == "C", "*"))

# filtering MT2 LTRs and MERVL-int
mesAur1_MT2_MERVL <- 
  mesAur1_LTRs %>%
  filter(grepl("MT2|MERVL", element_name))

# ploting width distribution
ggplot(data = filter(mesAur1_MT2_MERVL, grepl("MERVL", element_name)), aes(width)) + 
  geom_histogram(binwidth = 10)
