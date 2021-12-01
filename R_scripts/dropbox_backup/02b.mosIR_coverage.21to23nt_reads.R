### INFO: plot coverage over pCMV MosIR EGFP with 21-23 nt reads (normalized)
### DATE: Fri Jul 20 17:34:24 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working directory
setwd("/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_spleen.MosIR.small_RNAseq.2021_Jan/Analysis/MosIR_coverage")

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

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get coverage from bam file
coverageMosIR <- function(bam_path, which_gr, on_minus_strand){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = on_minus_strand))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
  bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
  
  # get length of chromosome on which is feature located
  seq_length <- seqlengths(bam_gr)[as.character(seqnames(which_gr))]
  
  # get coverage
  coverage_df <- 
    bam_gr %>% 
    coverage(.) %>% 
    .[unique(seqnames(bam_gr))] %>% 
    as(., "IntegerList") %>% 
    unlist(.) %>% 
    unname(.)
  
  if(length(coverage_df) == 0){
    coverage_df <- tibble(pos = 1:seq_length, 
                          coverage = 0)
  }else{
    coverage_df <- tibble(pos = 1:seq_length, 
                          coverage = coverage_df)
  }
  
  # set position to 0
  coverage_df %<>% 
    dplyr::filter((pos >= start(which_gr)) & (pos <= end(which_gr))) %>% 
    dplyr::mutate(pos = 1:nrow(.))
  
  # return 
  return(coverage_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# MosIR .bed path
bed_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/pCAG-EGFP_MosIR.coordinates.bed"

# experiment paths
experiment_path <- "/common/WORK/fhorvat/Projekti/Svoboda/siRNA.Valeria/datasets/mouse_spleen.MosIR.small_RNAseq.2021_Jan"

# mapped path
mapped_path <- file.path(experiment_path, "Data/Mapped/STAR_mm10")

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*\\.bam$", full.names = T)

# library size path
library_size_path <- file.path(mapped_path, "4_library_size", "library_sizes.txt")

######################################################## READ DATA
# read bed with coordinates
mosir_gr <- rtracklayer::import.bed(con = bed_path)
mosir_gr <- mosir_gr[mosir_gr$name != "EGFP"]

# read library size df
library_size_df <-  readr::read_delim(file = library_size_path, delim = "\t", col_names = c("sample_id", "library_size"))

######################################################## MAIN CODE
### prepare files
# set whole plasmid coordinates
plasmid_gr <- GenomicRanges::GRanges(seqnames = "pCAG-EGFP_MosIR",
                                     ranges = IRanges(start = 1, end = 6600),
                                     gene_id = "pCAG-EGFP_MosIR")

# clean library size
library_size_tidy <- 
  library_size_df %>% 
  dplyr::filter(str_detect(sample_id, "\\.19to32nt$")) %>% 
  dplyr::mutate(sample_id = str_remove(sample_id, "\\.19to32nt$"))


### create table with coverage
# loop through samples in experiment 
coverage_tb_full <- purrr::map(bam_paths, function(path){
  
  # get bam name
  bam_name <- basename(path) %>% str_remove(., ".bam")
  
  # print bam name
  cat(bam_name, "\n")
  
  # get coverage of both mosIR arms in one data.frame, normalize for library size
  coverage_tb <- 
    rbind(coverageMosIR(bam_path = path, which_gr = plasmid_gr, on_minus_strand = F) %>% 
            dplyr::mutate(strand = "plus"), 
          coverageMosIR(bam_path = path, which_gr = plasmid_gr, on_minus_strand = T) %>% 
            dplyr::mutate(strand = "minus", 
                          coverage = - coverage)) %>% 
    dplyr::mutate(sample_id = bam_name) %>%
    dplyr::left_join(., library_size_df, by = "sample_id") %>%
    dplyr::mutate(library_size = round((library_size / 1e6), 4),
                  rpm = coverage / library_size)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)

# create table for plot
plot_tb <- 
  coverage_tb_full %>% 
  dplyr::filter(str_detect(sample_id, "\\.SE$")) %>% 
  dplyr::mutate(genotype_dicer = str_extract(sample_id, "DicerSom|DicerWT|DicerX|WT"), 
                genotype_mosir = str_extract(sample_id, "MosIR"), 
                genotype_pkr = str_extract(sample_id, "delPKR")) %>% 
  dplyr::mutate(genotype_mosir = replace(genotype_mosir, is.na(genotype_mosir), "no_transfection"), 
                genotype_pkr = replace(genotype_pkr, genotype_dicer == "WT", "no_pkr"), 
                genotype_pkr = replace(genotype_pkr, is.na(genotype_pkr), "PKR")) %>% 
  dplyr::mutate(genotype_dicer = factor(genotype_dicer, levels = c("DicerSom", "DicerX", "DicerWT", "WT")), 
                genotype_mosir = factor(genotype_mosir, levels = c("MosIR", "no_transfection")), 
                genotype_pkr = factor(genotype_pkr, levels = c("PKR", "delPKR", "no_pkr"))) %>% 
  dplyr::arrange(genotype_dicer, genotype_mosir, genotype_pkr)


### plot whole plasmid
# plot
coverage_plot <-
  ggplot() +
  geom_rect(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(vars(sample_id), nrow = length(unique(plot_tb$sample_id))) + 
  scale_y_continuous(limits = c(-1.6, 1.6)) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(plot = coverage_plot,
       filename = file.path(outpath, str_c("MosIR_coverage.21to23nt_reads.whole_plasmid.RPM.png", sep = ".")),
       width = 8,
       height = 15)


### plot only MosIR repeat
# plot
coverage_plot <-
  ggplot() +
  geom_rect(data = plot_tb, aes(xmin = pos, xmax = pos + 1, ymin = 0, ymax = rpm, fill = strand)) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(vars(sample_id), nrow = length(unique(plot_tb$sample_id))) + 
  scale_x_continuous(limits = mosir_gr %>% asdf(.) %$% c(start, end) %>% range) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# save
ggsave(plot = coverage_plot,
       filename = file.path(outpath, str_c("MosIR_coverage.21to23nt_reads.repeat_only.RPM.png", sep = ".")),
       width = 8,
       height = 15)
