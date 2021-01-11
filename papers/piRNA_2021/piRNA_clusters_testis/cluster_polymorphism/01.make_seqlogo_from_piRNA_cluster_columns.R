### INFO: 
### DATE: Tue May 05 19:14:57 2020
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

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
library(rtracklayer)
library(GenomicAlignments)
library(Biostrings)
library(ggseqlogo)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# calculate proportion
proportion <- function(x){
  rs <- sum(x)
  return(x / rs)
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# small RNA-seq path
bams_path <- list.files(inpath, "\\.bam$", full.names = T)

######################################################## READ DATA
# read bams to list 
bams_list <- purrr::map(bams_path, function(path){
  
  # read bam, get sequences
  bam_gr <- GenomicAlignments::readGAlignments(file = path, use.names = TRUE, param = ScanBamParam(what = "seq"))
  
  # return
  return(bam_gr)
  
}) %>% 
  set_names(., basename(bams_path) %>% str_remove(., "\\.24to31nt\\.bam"))


######################################################## MAIN CODE
# combine 13 WT and 13 KO samples
bams_13dpp <- do.call(c, unname(bams_list[str_detect(names(bams_list), "13")]))

# get the stacked sequences
bams_stacked <- 
  bams_13dpp %>% 
  as_tibble(.) %>%
  dplyr::group_by(seqnames, start, end, strand) %>%
  dplyr::summarize(count = n()) %>% 
  arrange(desc(count)) %>% 
  dplyr::filter(count > 10) %>%
  GRanges(.)

# set names
names(bams_stacked) <- str_c(seqnames(bams_stacked), ":", start(bams_stacked), "-", end(bams_stacked), ":", strand(bams_stacked), " (", mcols(bams_stacked)$count, ")")


### get nucletoide frequency matrix for each stack
# create list
nucl_freq_list <- purrr::map(names(bams_stacked), function(name){
  
  # get reads in one stack
  bam_filt <- bams_13dpp[start(bams_13dpp) == start(bams_stacked[name]) & 
                           end(bams_13dpp) == end(bams_stacked[name]) &
                           strand(bams_13dpp) == strand(bams_stacked[name])]
  
  # get sequences cut based on CIGAR, calculate frequency
  nucl_freq <- 
    GenomicAlignments::sequenceLayer(x = mcols(bam_filt)$seq, cigar = cigar(bam_filt)) %>%
    Biostrings::consensusMatrix(., as.prob = T) %>% 
    .[rowSums(.) > 0, ]
  
  # return
  return(nucl_freq)
  
}) %>% 
  set_names(names(bams_stacked))

# plot
seqlogo_plot <- 
  ggplot() +
  geom_logo(nucl_freq_list, method = "custom") + 
  theme_logo() + 
  facet_wrap(~seq_group, ncol = 1, scales = "fixed")

# save 
ggsave(filename = file.path(outpath, str_c("testis_smallRNAseq.piRNA_cluster.13dpp_WT_KO.all.seqlogo.pdf")), 
       plot = seqlogo_plot, 
       width = 20, 
       height = 2.5 * length(nucl_freq_list), 
       limitsize = FALSE)


### get and save those which have polymorphism
# get names of all logos which have any polymorphism
stack_coords <- purrr::map(nucl_freq_list, function(tb) any(colSums(tb == 0) < 3)) %>% unlist %>% names(nucl_freq_list)[.]

# filter list
nucl_freq_list_polymorphisms <- nucl_freq_list[names(nucl_freq_list) %in% stack_coords]

# plot
seqlogo_plot <- 
  ggplot() +
  geom_logo(nucl_freq_list_polymorphisms, method = "custom") + 
  theme_logo() + 
  facet_wrap(~seq_group, ncol = 1, scales = "fixed")

# save 
ggsave(filename = file.path(outpath, str_c("testis_smallRNAseq.piRNA_cluster.13dpp_WT_KO.poly.seqlogo.pdf")), 
       plot = seqlogo_plot, 
       width = 20, 
       height = 2.5 * length(nucl_freq_list_polymorphisms), 
       limitsize = FALSE)

  
### get one stack and plot it in all samples
# get names of all logos which have any polymorphism
stack_coords <- purrr::map(nucl_freq_list, function(tb) any(colSums(tb == 0) < 3)) %>% unlist %>% names(nucl_freq_list)[.]

# plot them in 13dpp and 21 dpp WT sample
seqlogo_plots_list <- purrr::map(stack_coords, function(coords){
  
  # get stack coordinates
  stack_coords_tb <- 
    tibble(coords = coords) %>% 
    dplyr::mutate(coords = str_remove(coords, " .*")) %>% 
    tidyr::separate(coords, into = c("seqnames", "start", "end", "strand"), sep = "-|:")
  
  # get nucleotide freq in all bams
  nucl_freq_list_samples <- purrr::map(names(bams_list)[str_detect(names(bams_list), "13dpp|WT_21dpp")], function(name){
    
    # get one sample
    bam_sample <- bams_list[[name]]
    
    # get reads in one stack
    bam_filt <- bam_sample[start(bam_sample) == stack_coords_tb$start & 
                             end(bam_sample) == stack_coords_tb$end &
                             strand(bam_sample) == stack_coords_tb$strand]
    
    # get sequences cut based on CIGAR, calculate frequency
    if(length(bam_filt) > 0){
      
      # calculate frequency
      nucl_freq <- 
        GenomicAlignments::sequenceLayer(x = mcols(bam_filt)$seq, cigar = cigar(bam_filt)) %>%
        Biostrings::consensusMatrix(., as.prob = T) %>% 
        .[rownames(.) %in% c("A", "C", "G", "T"), ]
      
    }else{
      
      # create empty matrix
      nucl_freq <- matrix(rep(0.0001*4), nrow = 4, byrow = T)
      rownames(nucl_freq) <- c("A", "C", "G", "T")
      
    }
    
    # return
    return(list(nucl_freq = nucl_freq, length = length(bam_filt)))
    
  }) %>% 
    set_names(names(bams_list)[str_detect(names(bams_list), "13dpp|WT_21dpp")])
  
  
  # get lengths
  nucl_freq_list_samples_length <- purrr::map(nucl_freq_list_samples, 2)
  
  # get matrix
  nucl_freq_list_samples <- purrr::map(nucl_freq_list_samples, 1)
  
  # add length to names
  if(all(names(nucl_freq_list_samples) == names(nucl_freq_list_samples_length))){
    
    # paste 
    names(nucl_freq_list_samples) <- str_c(names(nucl_freq_list_samples), " (", unlist(nucl_freq_list_samples_length), ")")
    
  }else{
    
    # error
    stop("Something's wrong")
  }
  
  
  # add 13 dpp logo to list
  nucl_freq_list_samples <- c(nucl_freq_list[names(nucl_freq_list) == coords], nucl_freq_list_samples)
  names(nucl_freq_list_samples)[1] <- str_c("13 dpp WT and KO samples ", names(nucl_freq_list_samples)[1])
  
  # plot
  seqlogo_plot <- 
    ggplot() +
    geom_logo(nucl_freq_list_samples, method = "custom") + 
    theme_logo() + 
    facet_wrap(~seq_group, ncol = 1, scales = "fixed")
  
  # # save 
  # ggsave(filename = file.path(outpath, str_c("testis_smallRNAseq.piRNA_cluster.", 
  #                                            str_replace_all(coords, ":|-", "_") %>% str_remove(., "_\\+$|_-$"),
  #                                            ".seqlogo.pdf")), 
  #        plot = seqlogo_plot, 
  #        width = 15, 
  #        height = 2.5 * length(nucl_freq_list_samples), 
  #        limitsize = FALSE)
  
  # return
  return(seqlogo_plot)
  
})

# save
pdf(file = file.path(outpath, str_c("testis_smallRNAseq.piRNA_clusters", "individual_stacks", "seqlogo.pdf", sep = ".")), width = 15, height = 2.5 * 7)
invisible(purrr::map(seqlogo_plots_list, print))
dev.off()
