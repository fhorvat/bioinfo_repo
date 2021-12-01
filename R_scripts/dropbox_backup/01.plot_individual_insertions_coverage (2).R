### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/methylation/coverage")

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

library(ggrepel)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- file.path(getwd(), "..")

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi.trimmed"

# bigWigs path
bw_path_list <- list.files(mapped_path, ".*\\.deduplicated\\.raw\\.bw$", full.names = T) 

# methylation data path
meth_path_list <- list.files(mapped_path, ".*\\.deduplicated\\.bismark\\.cov\\.gz", full.names = T)

# bed path
bed_path <- "../.."
bed_path <- list.files(bed_path, ".*\\.FLI_elements\\.bed", full.names = T)

######################################################## READ DATA
# read bed
fli_gr <- rtracklayer::import.bed(bed_path)

# read coverage from bigWig
element_coverage_list <- purrr::map(bw_path_list, function(path){
  
  # read bigWig file
  coverage <- rtracklayer::import(path)
  
  # add genotype 
  mcols(coverage)$genotype <- str_extract(basename(path), "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO")
  
  # return
  return(coverage)
  
})

# read methylation data
meth_tb <- purrr::map(meth_path_list, function(path){
  
  # read table
  readr::read_delim(path, delim = "\t", col_names = c("seqnames", "start", "end", "percentage_meth", "count_meth", "count_unmeth")) %>% 
    dplyr::mutate(sample_id = basename(path) %>% str_remove(., "_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz|\\.merged\\.bismark\\.cov\\.gz"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO"))) %>% 
  tidyr::unite(col = coordinates, seqnames, start, end, sep = " ", remove = F) %>% 
  dplyr::select(seqnames, start, end, count_meth, count_unmeth, genotype, coordinates)

######################################################## MAIN CODE
### prepare files
# create element table
fli_tb <- 
  fli_gr %>% 
  as_tibble(.)

# join coverage for each genotype
element_coverage_list <- do.call(c, element_coverage_list)
element_coverage_list <- split(element_coverage_list, mcols(element_coverage_list)$genotype)

# join methylation for each element
element_methylation_list <- split(meth_tb, meth_tb$genotype)
element_methylation_list <- purrr::map(element_methylation_list, GRanges)


### find coverage for each element
element_coverage_tb <- purrr::map(names(element_coverage_list), function(genotype){
  
  # get one genotype
  element_coverage_full <- element_coverage_list[[genotype]]
  
  # find overlaps between FLIs and coverage
  hits <- findOverlaps(fli_gr, element_coverage_full, ignore.strand = T)
  
  # extract all overalaping features from subject as list
  element_coverage <- extractList(element_coverage_full, as(hits, "List"))
  
  # intersect with element coordinates
  element_coverage <- pintersect(element_coverage, fli_gr, ignore.strand = T)
  
  # set names
  names(element_coverage) <- mcols(fli_gr)$name 
  
  # set to table
  coverage_tb <- purrr::map(names(element_coverage), function(rmsk_id){
    
    cat(rmsk_id, "\n")
    
    # get FLI coordinates
    fli_coords <- fli_gr[mcols(fli_gr)$name == rmsk_id]
    
    # extract coverage
    cov_gr <- element_coverage[[rmsk_id]]
    
    # check if there is any coverage
    if(length(cov_gr) > 0){
      
      # add start and end from FLI coordinates
      if(start(cov_gr)[1] > start(fli_coords)){
        
        # add start
        start_gr <- GenomicRanges::restrict(fli_coords, end = (start(cov_gr)[1] - 1))
        mcols(start_gr) <- mcols(cov_gr[1])
        mcols(start_gr)$score <- 0
        cov_gr <- c(start_gr, cov_gr)
        
      }
      
      if(end(cov_gr)[length(cov_gr)] < end(fli_coords)){
        
        # add end
        end_gr <- GenomicRanges::restrict(fli_coords, start = (end(cov_gr)[length(cov_gr)] + 1))
        mcols(end_gr) <- mcols(cov_gr[length(cov_gr)])
        mcols(end_gr)$score <- 0
        cov_gr <- c(cov_gr, end_gr)
        
      }
      
    }else{
      
      # create empty range
      cov_gr <- fli_coords
      mcols(cov_gr) <- NULL
      mcols(cov_gr)$score <- 0
      mcols(cov_gr)$genotype <- genotype
      mcols(cov_gr)$hit <- TRUE
      
    }
    
    # to table 
    cov_tb <- 
      cov_gr %>% 
      as_tibble(.) %>% 
      dplyr::mutate(strand = as.character(strand(fli_coords)), 
                    rmsk_id = rmsk_id) %>% 
      dplyr::select(-hit)
    
    # return
    return(cov_tb)
    
  }) %>% 
    dplyr::bind_rows(.)
  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO")))


### find methylation for each element 
element_methylation_tb <- purrr::map(names(element_methylation_list), function(genotype){
  
  # get one genotype
  element_methylation_full <- element_methylation_list[[genotype]]
  
  # find overlaps between FLIs and coverage
  hits <- findOverlaps(fli_gr, element_methylation_full, ignore.strand = T)
  
  # extract all overalaping features from subject as list
  element_methylation <- extractList(element_methylation_full, as(hits, "List"))
  
  # set names
  names(element_methylation) <- mcols(fli_gr)$name 
  
  # unlist, set names as rmsk_id
  element_methylation_tb <- 
    element_methylation %>% 
    unlist(.)
  mcols(element_methylation_tb)$rmsk_id <- names(element_methylation_tb)
  names(element_methylation_tb) <- NULL
  
  # to table
  element_methylation_tb <- 
    element_methylation_tb %>% 
    as_tibble(.) %>% 
    dplyr::select(-strand) %>% 
    dplyr::left_join(., fli_tb %>% dplyr::select(rmsk_id = name, strand), by = "rmsk_id")

  # return 
  return(element_methylation_tb)
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(meth_percentage = 100*(count_meth / (count_meth + count_unmeth)))


### plot
purrr::map(unique(fli_tb$name), function(fli_id){
  
  # filter tables
  fli_tb_filt <- fli_tb %>% dplyr::filter(name == fli_id)
  cov_plot_tb <- element_coverage_tb %>% dplyr::filter(rmsk_id == fli_id)
  meth_plot_tb <- element_methylation_tb %>% dplyr::filter(rmsk_id == fli_id)
  
  # coverage plot
  coverage_plot <- 
    ggplot() + 
    geom_rect(data = cov_plot_tb, aes(xmin = start, xmax = end, ymin = 0, ymax = score), fill = "black") +
    geom_rect(data = meth_plot_tb, aes(xmin = start, xmax = start + 1, ymin = 0, ymax = count_meth), color = "red", alpha = 1) +
    # geom_text_repel(data = meth_plot_tb, aes(x = start, y = count_meth, label = meth_percentage), 
    #                 force = 2, angle = 90, direction = "x", segment.color = NA) +
    facet_grid(rows = vars(genotype)) + 
    ggtitle(str_c(fli_tb_filt$name, fli_tb_filt$seqnames, fli_tb_filt$start, fli_tb_filt$end, fli_tb_filt$strand, sep = " ")) +
    ylab("") +
    xlab("") +
    # scale_y_continuous(limits = c(0, 90)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          # plot.title = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # reverse if on minus strand
  if(fli_tb_filt$strand == "-"){
    
    coverage_plot <- 
      coverage_plot +  
      scale_x_reverse(limits = c(fli_tb_filt$end, fli_tb_filt$start))
    
  } else{
    
    coverage_plot <- 
      coverage_plot +  
      scale_x_continuous(limits = c(fli_tb_filt$start, fli_tb_filt$end))
    
  }
  
  # save plot
  ggsave(plot = coverage_plot,
         filename = file.path(outpath, str_c(str_c(fli_tb_filt$seqnames, fli_tb_filt$start, fli_tb_filt$end, sep = " "), ".jpg")),
         width = 15, height = 7.5)
  
  # return
  return(rmsk_id)
  
})
