### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
# options(bitmapType = "cairo")
wideScreen()
# Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1Ijoid2lja2VybWFuOTk5IiwiYSI6ImNqbjB4ZmdxdjFleWgzdnFvanl6ZHowbGIifQ.giaFevJaFe7KJ5_5X0j_PA')

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/phasing/NIH3T3_transfected.2018")

######################################################## LIBRARIES
library(dplyr)
library(readr)
library(magrittr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)

library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(Biostrings)
library(plotly)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# get coverage from bam file
phasingRegister <- function(bam_path, which_gr){
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = which_gr, 
                                                                flag = scanBamFlag(isMinusStrand = T))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
  bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
  
  # get starting positions of reads, normalize to start of mosIR arm, get mod 21
  if(as.character(strand(which_gr)) == "+"){
    
    # if read is on MosIR arm 1 just set start of arm 1 to 1
    values(bam_gr)$start_norm <- start(bam_gr) - start(which_gr) + 1
    values(bam_gr)$end_norm <- end(bam_gr) - start(which_gr) + 1
    
  }else{
    if(as.character(strand(which_gr)) == "-"){
      
      # if read is on MosIR arm 2 then start is really an end of mosIR arm 2
      values(bam_gr)$start_norm <- end(which_gr) - end(bam_gr) + 1
      values(bam_gr)$end_norm <- end(which_gr) - start(bam_gr) + 1
      
    }else{
      
      # if non true then input is wrong
      stop("Something went wrong!")
      
    }
  }
  
  ### get mod 22 of read starts
  read_start_mod <- (values(bam_gr)$start_norm %% 22)
  
  # table for plot
  read_register_df <- 
    tibble(read_register = read_start_mod) %>% 
    dplyr::group_by(read_register) %>% 
    dplyr::summarize(count = n()) %>% 
    dplyr::mutate(percent = (count / sum(count)) * 100) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-count) %>% 
    dplyr::right_join(tibble(read_register = 0:21), by = "read_register") %>% 
    dplyr::mutate(percent = replace(percent, is.na(percent), 0))
  
  # return 
  return(read_register_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# mapped path
mapped_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" 

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)

# get path of inverted repeat coordinates from mRNA
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/mm10.all_plasmids/IR_mRNAs.plasmid_coordinates/IR_mRNA.plasmid_coordinates.clean.csv"

######################################################## READ DATA
# read coordinates of inverted repeats from mRNA
ir_gr <- readr::read_csv(file = ir_path)

######################################################## MAIN CODE
# prepare mRNA IR coordinates on plasmids
ir_gr %<>% GenomicRanges::GRanges(.)
mcols(ir_gr) <- mcols(ir_gr)[, c("arm")]
names(mcols(ir_gr)) <- "gene_id"
mcols(ir_gr)$gene_biotype <- as.character(seqnames(ir_gr))
ir_gr <- GenomicRanges::reduce(ir_gr, min.gapwidth = 1000, ignore.strand = T)
strand(ir_gr) <- "-"

# loop through sample in experiment 
plot_facet_df <- purrr::map(.x = bam_paths, function(bam_path){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam"); cat(bam_name, "\n")

  # get IR name
  ir_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")
  
  # subset IR GenomicRanges
  ir_filt <- ir_gr[seqnames(ir_gr) == ir_name]
  
  # get coverage
  plot_df <-
    phasingRegister(bam_path = bam_path, which_gr = ir_filt) %>%
    dplyr::mutate(arm = as.character(seqnames(ir_filt)), 
                  sample_id = bam_name) %>% 
    dplyr::select(sample_id, arm, everything())
  
}) %>% 
  dplyr::bind_rows(.)

### plot as separate samples
# split data.frame to list
plot_list <- 
  plot_facet_df %>% 
  split(., plot_facet_df$sample_id)

# plot in a loop
purrr::map(names(plot_list), function(bam_name){
  
  bam_name <- "s_pU6_Lin28IR_r2.SE.mis_0"
  
  # print bam name
  cat(bam_name, "\n")
  
  # get plot data.frame
  plot_df <- 
    plot_list[[bam_name]] %>% 
    dplyr::mutate(read_register = replace(read_register, read_register == 0, 22)) %>% 
    dplyr::arrange(read_register) %>% 
    dplyr::bind_rows(., .[1, ]) %>% 
    dplyr::mutate(read_register = str_c("r", read_register))

  p <- 
    plot_ly(plot_df, 
            type = "scatterpolar", 
            mode = "lines",
            r = ~percent, 
            theta = ~read_register, 
            fill = "toself", 
            color = ~arm, 
            colors = "#F8766D") %>% 
    layout(polar = list(angularaxis = list(rotation = 90, direction = "clockwise", 
                                           tickmode = "array", 
                                           tickvals = plot_df$read_register, 
                                           ticktext = str_remove(plot_df$read_register, "r"), 
                                           tickfont = list(size = 19)), 
                        radialaxis = list(visible = T, range = c(0, 50),
                                          tickmode = "array", 
                                          tickvals = c(0, 25, 50),
                                          tickprefix = "     ", 
                                          tickfont = list(size = 19))))
  
  export(p, file.path(outpath, "one_strand", str_c(bam_name, "IR_phase.png", sep = ".")))
  
  # return
  return(bam_name)

})
