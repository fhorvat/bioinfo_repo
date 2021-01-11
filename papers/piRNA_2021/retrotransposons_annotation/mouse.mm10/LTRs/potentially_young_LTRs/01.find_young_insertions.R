### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/mouse.mm10/LTRs/potentially_young_LTRs")

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
library(GenomicFeatures)
library(rtracklayer)
library(DESeq2)

library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(stringdist)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# set ensembl version
ensembl_version <- 99

# genome path
genome_dir <- "/common/DB/genome_reference/mouse/mm10.GRCm38.GCA_000001635.2"

# joined repeatMasker path
rmsk_path <- file.path(genome_dir, "rmsk.mm10.20180919.joined_rmsk_id.fa.out.gz")

# clean repeatMasker path
rmsk_clean_path <- file.path(genome_dir, "rmsk.mm10.20180919.clean.fa.out.gz")

######################################################## READ DATA
# read joined repeatMasker
rmsk_tb <- readr::read_delim(rmsk_path, delim = "\t")

# read clean repeatMasker
rmsk_clean <- 
  readr::read_delim(rmsk_clean_path, delim = "\t") %>% 
  dplyr::rename(rmsk_id = rmsk_ID)

######################################################## MAIN CODE
# get widths
rmsk_tb %<>%    
  # dplyr::mutate(strand = ifelse(!(strand %in% c("+", "-")), "+", strand)) %>%
  # GRanges(.) %>%
  # as_tibble(.) %>%
  dplyr::mutate(rmsk_id = as.character(rmsk_id), 
                width = end - start + 1)

# filter LTRs
ltr_tb <- 
  rmsk_tb %>%
  dplyr::filter(repClass == "LTR") %>% 
  dplyr::filter(str_detect(repName, ".*/.*/")) %>% 
  dplyr::filter(insertion_class %in% c("whole", "within")) %>% 
  dplyr::mutate(first_ltr = repName %>% str_remove(., "/.*"), 
                last_ltr = repName %>% str_remove(., ".*/")) %>% 
  dplyr::filter(first_ltr == last_ltr, 
                !str_detect(first_ltr, "int"), 
                !str_detect(last_ltr, "int")) %>% 
  dplyr::select(seqnames, start, end, strand, width, repName, ltr_name = first_ltr, repFamily, insertion_class, rmsk_id) %>% 
  dplyr::mutate(ltr_class = str_extract(ltr_name, "^ERV|^IAP|^MLT|^MMERVK|^MT|^ORR|^RMER|^RLTR|^MER|^LTR|Gyp|^MRLTR|^BGLII|^MMERGLN|^MURVY")) %>% 
  dplyr::filter(width < 10000)

# # plot histogram of widths
# hist_width <-
#   ggplot(data = ltr_tb, aes(width, color = ltr_name)) +
#   geom_histogram(binwidth = 4) +
#   facet_wrap(vars(ltr_class), nrow = 6, ncol = 2) +
#   # geom_freqpoly(data = iap_tb, aes(width), binwidth = 10, size = 1.5) +
#   # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme_bw(base_size = 10) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.key = element_blank(),
#         plot.title = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   theme(legend.position = "none")
# 
# # save
# ggsave(filename = file.path(outpath, "LTRs.width_hist.facet.png"), plot = hist_width, width = 10, height = 30)


### get individual LTRs from full insertions
# get first and last LTR in insertion
ltrs_flanks_tb <- 
  rmsk_clean %>% 
  dplyr::filter(rmsk_id %in% ltr_tb$rmsk_id) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::filter(length(unique(strand)) == 1) %>% 
  dplyr::filter(start == max(start) | start == min(start)) %>% 
  dplyr::ungroup(.)

# get sequences of first and last LTR
ltrs_flanks_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, GRanges(ltrs_flanks_tb))

# add sequence to table, split to list
ltrs_tb <- 
  ltrs_flanks_tb %>%
  dplyr::mutate(ltr_seq = as.character(ltrs_flanks_seq), 
                rmsk_id = as.character(rmsk_id)) %>% 
  dplyr::group_by(rmsk_id) %>% 
  dplyr::arrange(start) %>% 
  dplyr::mutate(ltr_in_insert = ifelse(strand == "+", c("5pLTR", "3pLTR"), rev(c("5pLTR", "3pLTR")))) %>% 
  dplyr::ungroup(.) %>% 
  tidyr::unite(coordinates, seqnames, start, end, sep = " ") %>% 
  dplyr::select(rmsk_id, ltr_in_insert, coordinates, sequence = ltr_seq) %>% 
  tidyr::pivot_wider(id_cols = "rmsk_id", names_from = ltr_in_insert, values_from = c(coordinates, sequence)) %>% 
  dplyr::mutate(width_3pLTR = nchar(sequence_3pLTR), 
                width_5pLTR = nchar(sequence_5pLTR)) %>% 
  dplyr::mutate(align_score = pairwiseAlignment(sequence_5pLTR, sequence_3pLTR, scoreOnly = T), 
                edit_distance = stringdist::stringdist(sequence_5pLTR, sequence_3pLTR, method = "lv")) %>% 
  dplyr::left_join(., ltr_tb %>% 
                     tidyr::unite(insertion_coordinates, seqnames, start, end, sep = " ") %>% 
                     dplyr::select(insertion_coordinates, strand, repName, ltr_class, ltr_name, rmsk_id, insertion_class), by = "rmsk_id") %>% 
  dplyr::select(insertion_coordinates, strand, repName, ltr_class, ltr_name, rmsk_id, insertion_class, align_score, edit_distance, 
                coordinates_5pLTR, coordinates_3pLTR, width_5pLTR, width_3pLTR, sequence_5pLTR, sequence_3pLTR)

# filter
ltrs_tb_filt <- 
  ltrs_tb %>% 
  # dplyr::filter(align_score > 0) %>% 
  dplyr::filter(edit_distance <= 10) %>% 
  dplyr::filter(abs(width_5pLTR - width_3pLTR) <= 10) %>% 
  dplyr::arrange(ltr_class, repName)


### plot width distribution
# get coordinates from table
ltrs_coords <-
  ltrs_tb_filt %>%
  dplyr::select(rmsk_id, ltr_class, ltr_name, coordinates_5pLTR, coordinates_3pLTR) %>%
  tidyr::separate(coordinates_5pLTR, into = c("seqnames_5pLTR", "start_5pLTR", "end_5pLTR"), sep = " ") %>%
  tidyr::separate(coordinates_3pLTR, into = c("seqnames_3pLTR", "start_3pLTR", "end_3pLTR"), sep = " ") %>%
  dplyr::mutate(start_5pLTR = as.numeric(start_5pLTR), end_5pLTR = as.numeric(end_5pLTR),
                start_3pLTR = as.numeric(start_3pLTR), end_3pLTR = as.numeric(end_3pLTR)) %>%
  dplyr::select(-c(seqnames_5pLTR, seqnames_3pLTR)) %>%
  dplyr::mutate(ltr_distance = abs(end_5pLTR - start_3pLTR))

# plot histogram of LTR distances
hist_width <-
  ggplot(data = ltrs_coords, aes(ltr_distance, color = ltr_class)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(vars(ltr_class), nrow = 7) +
  # geom_freqpoly(data = iap_tb, aes(width), binwidth = 10, size = 1.5) +
  # scale_x_continuous(limits = c(100, 10000), breaks = seq(0, 10000, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

# save
ggsave(filename = file.path(outpath, "LTRs.distances_hist.png"), plot = hist_width, width = 10, height = 10)


# ### plot individual LTR classes
# # split by LTR class
# ltrs_coords_list <- split(ltrs_coords, ltrs_coords$ltr_class)
# 
# # plot
# purrr::map(ltrs_coords_list, function(ltr_tb){
#   
#   # plot histogram of LTR distances
#   hist_width <-
#     ggplot(data = ltr_tb, aes(ltr_distance, fill = ltr_name, color = ltr_name)) +
#     geom_histogram(binwidth = 1) +
#     # geom_freqpoly(data = ltr_tb, aes(ltr_distance, color = ltr_name), binwidth = 5, size = 1.5) +
#     scale_x_continuous() +
#     expand_limits(x = c(0, max(ltr_tb$ltr_distance) + 10)) + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     theme_bw(base_size = 10) +
#     theme(legend.position = "bottom",
#           legend.title = element_blank(),
#           legend.background = element_blank(),
#           legend.key = element_blank(),
#           plot.title = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()) +
#     theme(legend.position = "bottom")
#   
#   # save
#   ggsave(filename = file.path(outpath, str_c("LTRs.distances_hist.", unique(ltr_tb$ltr_class), ".png")), plot = hist_width, width = 10, height = 10)
#   
#   # return 
#   return(unique(ltr_tb$ltr_class))
#   
# })


### save table and fasta
# save 
readr::write_csv(ltrs_tb_filt, path = file.path(outpath, str_c("potentially_young_LTRs.edit_distance", 10, "csv", sep = ".")))

# get sequences
ltrs_seq <- 
  ltrs_tb_filt %>% 
  tidyr::separate(col = insertion_coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  GRanges(.) %>% 
  getSeq(BSgenome.Mmusculus.UCSC.mm10, .)

# add to table, split
ltrs_tb_split <- 
  ltrs_tb_filt %>% 
  dplyr::mutate(full_seq = as.character(ltrs_seq)) %>% 
  split(., .$ltr_class)

# transform to DNAStringSet, save as fasta
purrr::map(ltrs_tb_split, function(ltr){
  
  # get DNAStringSet
  ltr_seq <- Biostrings::DNAStringSet(ltr$full_seq)
  names(ltr_seq) <- str_c(ltr$repName, ltr$rmsk_id, sep = ".")
  
  # save 
  Biostrings::writeXStringSet(x = ltr_seq, 
                              filepath = file.path(outpath, str_c("potentially_young_LTRs.edit_distance", 10, "full_insertion_sequences", unique(ltr$ltr_class), "fasta", sep = ".")))
  
  # return 
  return(unique(ltr$ltr_class))
  
})

