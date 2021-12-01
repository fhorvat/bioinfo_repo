### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/methylation/merge_CpG_based_on_coverage")

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
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/hamster_KO/datasets/hamster_oocyte_Mov10l.bisulfite/Data/Mapped/Bismark_Siomi"

# bigWigs path
bw_path_list <- list.files(mapped_path, ".*\\.deduplicated\\.raw\\.bw$", full.names = T) 

# methylation data path
meth_path_list <- list.files(mapped_path, ".*\\.deduplicated\\.bismark\\.cov\\.gz", full.names = T)

######################################################## READ DATA
# read coverage from bigWig
coverage_list <- purrr::map(bw_path_list, function(path){
  
  # read bigWig file
  coverage <- rtracklayer::import(path)
  
  # add genotype 
  mcols(coverage)$genotype <- str_extract(basename(path), "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO")
  
  # return
  return(coverage)
  
})

# read methylation data
meth_gr <- purrr::map(meth_path_list, function(path){
  
  # read table
  readr::read_delim(path, delim = "\t", col_names = c("seqnames", "start", "end", "percentage_meth", "count_meth", "count_unmeth")) %>% 
    dplyr::mutate(sample_id = basename(path) %>% 
                    str_remove(., "_bismark_bt2_pe\\.deduplicated\\.bismark\\.cov\\.gz|\\.merged\\.bismark\\.cov\\.gz"))
  
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::mutate(genotype = str_extract(sample_id, "Mov10l1_WT|Mov10l1_HET|Mov10l1_KO"),
                genotype = factor(genotype, levels = c("Mov10l1_WT", "Mov10l1_HET", "Mov10l1_KO"))) %>% 
  tidyr::unite(col = coordinates, seqnames, start, end, sep = " ", remove = F) %>% 
  dplyr::select(seqnames, start, end, coordinates, count_meth, count_unmeth, genotype) %>% 
  GRanges(.)

######################################################## MAIN CODE
### clean and prepare files
# join coverage for each genotype
coverage_list <- do.call(c, coverage_list)
coverage_list <- split(coverage_list, mcols(coverage_list)$genotype)

# join methylation for each genotype
meth_list <- split(meth_gr, mcols(meth_gr)$genotype)


### merge coverage >= 4 reads
# separately for each genome
element_coverage_tb <- purrr::map(names(element_coverage_list), function(genotype){
  
  ### clean data
  # get methylation sites in genotype
  meth_gr <- meth_list[[genotype]]
  
  # get CpGs supported by at least four reads, join together
  raw_coverage <- 
    coverage_list[[genotype]] %>% 
    .[mcols(.)$score >= 4] %>% 
    GenomicRanges::reduce(., ignore.strand = T) 
  
  
  ## classify each CpG site
  # find overlaps between CpGs and coverage
  hits <- findOverlaps(meth_gr, raw_coverage, ignore.strand = T)
  
  # extract merged coverage
  coverage_hits <- 
    raw_coverage[subjectHits(hits)] %>% 
    as_tibble(.) %>% 
    tidyr::unite(coverage_coordinates, seqnames, start, end, sep = " ") %>% 
    dplyr::select(coverage_coordinates)
  
  # extract all CpGs, summarize by coverage hits
  meth_hits <- 
    meth_gr[queryHits(hits)] %>% 
    as_tibble(.) %>% 
    dplyr::mutate(coverage_coordinates = coverage_hits$coverage_coordinates) %>% 
    dplyr::group_by(coverage_coordinates) %>% 
    dplyr::summarise(count_meth = sum(count_meth),
                     count_unmeth = sum(count_unmeth)) %>% 
    dplyr::ungroup(.) %>% 
    dplyr::mutate(genotype = genotype) %>% 
    tidyr::separate(coverage_coordinates, c("seqnames", "start", "end"), sep = " ", remove = F) %>% 
    dplyr::select(seqnames, start, end, coverage_coordinates, count_meth, count_unmeth, genotype) %>% 
    GRanges(.)
    
})





# set to table
coverage_tb <- purrr::map(names(element_coverage), function(rmsk_id){
  
  rmsk_id <- "190629"
  
  # get FLI coordinates
  fli_coords <- fli_gr[mcols(fli_gr)$name == rmsk_id]
  
  # extract coverage
  cov_gr <- element_coverage[[rmsk_id]]
  
  # add start and end
  if(start(cov_gr)[1] > start(fli_coords)){
    
    start_gr <- GenomicRanges::restrict(fli_coords, end = (start(cov_gr)[1] - 1))
    mcols(start_gr) <- mcols(cov_gr[1])
    mcols(start_gr)$score <- 0
    
  }
  
  as_tibble(.) %>% 
    dplyr::mutate(start_in_element = start - start(fli_coords), 
                  end_in_element = end - start(fli_coords))
  
  as(., "IntegerList") %>% 
    unlist(.)
  tibble(coverage = .) %>% 
    dplyr::mutate(repName = seqname,
                  sample_id = basename(path), 
                  pos = 1:nrow(.)) %>% 
    dplyr::select(repName, sample_id, pos, coverage)
  
}) %>% 
  dplyr::bind_rows(.)


# read coverage data
coverage_tb_full <- purrr::map(bw_path, function(path){
  
  # read bw
  bw_gr <- rtracklayer::import.bw(path, as = "GRanges")
  bw_gr <- subsetByOverlaps(bw_gr, fli_gr, ignore.strand = T)
  
  # find coverage of individual FLI
  hits <- findOverlaps(fli_gr, bw_gr, ignore.strand = T)
  
  # extract all overalaping features from subject as list
  element_coverage <- extractList(bw_gr, as(hits, "List"))
  
  # get coverage
  coverage_list <- 
    bw_gr %>% 
    as(., "IntegerList") %>% 
    unlist(.) %>% 
    split(., names(.)) 
  

  
  # return
  return(coverage_tb)
  
}) %>% 
  dplyr::bind_rows(.)

### get data for plot 
# filter table
coverage_tb <- 
  coverage_tb_full %>% 
  dplyr::mutate(sample_id = str_remove_all(sample_id, "\\.converted_bismark_bt2.raw_coverage.bw"))

# filter methylation table - don't take in account CpG sites with less than 30 reads
meth_tb <- 
  meth_tb_full %>% 
  dplyr::filter((methylated + unmethylated) >= 30)

# polygon total coverage
pg <- tibble(px = coverage_tb$pos, 
             py = coverage_tb$coverage, 
             sample_id = coverage_tb$sample_id, 
             repName = coverage_tb$repName) 

# polygon methylation
pg_meth <- tibble(px = meth_tb$start, 
                  percentage = meth_tb$percentage, 
                  sample_id = meth_tb$sample_id, 
                  repName = meth_tb$repName)

# loop through individual insertions
for(insertion in unique(pg$repName)){
  
  # filter coverage table
  pg_filt <- 
    pg %>% 
    dplyr::filter(repName == insertion)
  
  # filter methylation table
  pg_meth_tb <- 
    pg_meth %>% 
    dplyr::filter(repName == insertion) %>% 
    dplyr::mutate(py = (percentage * 0.01 * max(pg_filt$py)), 
                  percentage = round(percentage, 2))
  
  # create plot
  baseplot <- 
    ggplot() +
    geom_rect(data = pg_filt, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), fill = "black") +
    geom_rect(data = pg_meth_tb, aes(xmin = px, xmax = px + 1, ymin = 0, ymax = py), color = "red", alpha = 1) +
    geom_text_repel(data = pg_meth_tb, aes(x = px, y = py, label = percentage), force = 2, angle = 90, direction = "x", segment.color = NA) +
    facet_grid(cols = vars(sample_id), 
               scales = "free_y") +
    ggtitle(str_replace_all(insertion, ":|-", " ")) + 
    ylab("") +
    xlab("") +
    # scale_x_continuous(limits = c(0, 6500)) +
    # scale_y_continuous(limits = c(0, 90)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          # plot.title = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # save plot
  ggsave(plot = baseplot,
         filename = file.path(outpath, "test.png"),
         width = 15, height = 7.5)
  
  # # save plot
  # ggsave(plot = baseplot, 
  #        filename = file.path(outpath, str_c("raw_coverage", insertion, "bw", "jpg", sep = ".")), 
  #        width = 8, height = 2)
  
}
