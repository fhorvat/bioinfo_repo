### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2020_Apr.brain_spleen/Analysis/virus_genome_coverage")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(data.table)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

library(ggridges)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get coverage from bam file
extractCigar <- function(bam_path, which_gr){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., "\\.bam")
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr)) %>% 
    unlist(.)

  # create cigar table, summarize
  cigar_tb <-
    tibble(cigar = cigar(bam_gr)) %>% 
    dplyr::count(cigar) %>% 
    dplyr::mutate(sample_id = bam_name) %>% 
    dplyr::select(sample_id, cigar, count = n)
  
  # return
  return(cigar_tb)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/RNAi.Eliska/DicerX_viral_infection/datasets/2020_Apr.brain_spleen/Data/Mapped/STAR_mm10"

######################################################## READ DATA

######################################################## MAIN CODE
# set whole plasmid coordinates
viral_gr <- GenomicRanges::GRanges(seqnames = "NC_001672.1",
                                   ranges = IRanges(start = 1, end = 11141),
                                   gene_id = "NC_001672.1")

### get viral genome counts for all samples, normalize, plot
# list bams 
bam_paths <- list.files(path = mapped_path, pattern = "s_brain.*\\.SE\\.bam$", full.names = T)

# get library sizes 
library_tb <- 
  purrr::map(bam_paths, extractCigar, which_gr = viral_gr) %>% 
  dplyr::bind_rows(.)


### prepare data
# clean library
library_tidy <- 
  library_tb %>% 
  dplyr::filter(!str_detect(cigar, "I|D")) %>% 
  dplyr::mutate(alignment_length = str_extract(cigar, "[0-9]{2}M") %>% str_remove(., "M"), 
                left_soft = str_extract(cigar, "^[0-9]+S(?=[0-9]{2}M)") %>% str_remove(., "S"), 
                right_soft = str_extract(cigar, "(?<=[0-9]{2}M)[0-9]+S$") %>% str_remove(., "S"))

# get total library size
library_total <- 
  library_tidy %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(library_size = sum(count)) %>% 
  dplyr::ungroup(.) 


### get whole distribution
# get read length distribution
read_length_tb <- 
  library_tidy %>% 
  dplyr::select(sample_id, alignment_length, count) %>% 
  dplyr::mutate(alignment_length = as.numeric(alignment_length)) %>% 
  dplyr::group_by(sample_id, alignment_length) %>% 
  dplyr::summarise(count = sum(count)) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(count_total = sum(count)) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::mutate(percentage = round((count / count_total), 3) * 100)
# dplyr::filter(alignment_length <= 35)

# plot
read_length_plot <- 
  ggplot() +
  geom_histogram(data = read_length_tb %>% dplyr::filter(alignment_length <= 35), 
                 mapping = aes(x = alignment_length, y = percentage, fill = sample_id), width = 0.9, stat = "identity") +
  # scale_fill_viridis_d() +
  facet_grid(sample_id ~ .) +
  scale_x_discrete(limits = 16:35, breaks = 16:35, labels = as.character(16:35)) +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  theme(strip.text.y = element_text(size = 6, colour = "black"),
        # strip.background = element_blank(),
        strip.placement = "outside") +
  ggsave(filename = file.path(outpath, "read_length.viral_genome.histogram.limit.png"), width = 10, height = 15)


# Joy Division plot
ridge_plot <-
  ggplot() +
  geom_density_ridges(data = read_length_tb, mapping = aes(x = alignment_length, y = sample_id, height = percentage, fill = sample_id), scale = 1.5, stat = "identity") +
  scale_fill_viridis_d() +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(16:35, 58:80), labels = as.character(c(16:35, 58:80))) +
  theme_ridges() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ggsave(filename = file.path(outpath, "read_length.joy_division.perfect.png"), width = 15, height = 12)



### get summarized read lenghts
# get read length distribution
read_length_tb <- 
  library_tidy %>% 
  dplyr::select(sample_id, alignment_length, count) %>% 
  dplyr::mutate(alignment_length = as.numeric(alignment_length)) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(count_16to20nt = sum(count[alignment_length >=16 & alignment_length <= 20]),
                   count_21to23nt = sum(count[alignment_length >=21 & alignment_length <= 23]),
                   count_24to31nt = sum(count[alignment_length >=24 & alignment_length <= 31]),
                   count_32plus = sum(count[alignment_length > 31]),
                   count_total = sum(count)) %>%
  dplyr::ungroup(.) %>% 
  tidyr::pivot_longer(names_to = "length", values_to = "count", cols = -c(sample_id, count_total), names_prefix = "count_") %>% 
  dplyr::mutate(percentage = round((count / count_total), 3) * 100) %>% 
  dplyr::mutate(length = factor(length, levels = c("32plus", "24to31nt", "21to23nt", "16to20nt"), 
                                labels = c("32+", "24-31", "21-23", "16-20"))) %>% 
  dplyr::arrange(sample_id, desc(length)) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(pos = cumsum(percentage) - (0.5 * percentage)) %>%
  dplyr::ungroup(.) 

# plot
read_length_plot <- 
  ggplot() +
  geom_bar(data = read_length_tb, 
           mapping = aes(x = sample_id, y = percentage, fill = length), width = 0.8, stat = "identity") +
  geom_text(data = read_length_tb,
            mapping = aes(x = sample_id, y = pos, label = str_c(percentage, "%")), size = 4, stat = "identity") +
  geom_text(data = read_length_tb[read_length_tb$length == "32+", ],
            mapping = aes(x = sample_id, y = pos, label = str_c(percentage, "%")), size = 4, stat = "identity", color = "white") +
  scale_fill_viridis_d() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  ggsave(filename = file.path(outpath, "read_length.viral_genome.stacked_barplot.png"), width = 10, height = 20)


