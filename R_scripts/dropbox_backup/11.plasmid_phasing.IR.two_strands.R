### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
# options(bitmapType = "cairo")
wideScreen()

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
                                                                flag = scanBamFlag(isMinusStrand = NA))) %>% 
    unlist(.) %>% 
    .[str_detect(cigar(.), "21M|22M|23M")] %>% 
    as.data.frame(.) %>% 
    as.tibble(.)
  
  # separate reads on plus and minus strand
  read_register_df <- 
    bam_gr %>% 
    dplyr::mutate(start_norm = start - start(which_gr) + 1, 
                  end_norm = end - start(which_gr) + 1, 
                  read_register =  start_norm %% 22) %>% 
    dplyr::group_by(read_register, strand) %>% 
    dplyr::summarize(count = n()) %>% 
    dplyr::group_by(strand) %>% 
    dplyr::mutate(percent = (count / sum(count)) * 100)  %>% 
    dplyr::ungroup(.) %>% 
    dplyr::select(-count) %>% 
    dplyr::right_join(tibble(read_register = rep(0:21, each = 2), strand = rep(c("+", "-"), 22)), 
                      by = c("read_register", "strand")) %>%
    dplyr::mutate(percent = replace(percent, is.na(percent), 0))
  
  # return 
  return(read_register_df)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

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

# experiment paths
experiment_path <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.plasmids" %>% 
  magrittr::set_names(., str_remove(., "Data.*$") %>% basename(.))

# list of experiments
experiment_list <- names(experiment_path)

# loop through experiments
experiment <- experiment_list[1]

# print experiment
cat("Plotting experiment:", experiment, "\n")

# get mapped path 
mapped_path <- experiment_path[experiment]

# list bams 
bam_paths <- list.files(path = mapped_path, pattern = ".*mis_0.*bam$", full.names = T)

# loop through sample in experiment 
plot_facet_df <- purrr::map(.x = bam_paths, function(bam_path){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam"); cat(bam_name, "\n")
  
  # get IR name
  ir_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")
  
  # subset IR GenomicRanges
  ir_filt <- ir_gr[seqnames(ir_gr) == ir_name]
  ir_filt <- ir_filt[str_detect(mcols(ir_filt)$gene_id, "antisense_arm")]
  
  # get phasing of + and - strand of antisense arm 
  phasing_df <- 
    phasingRegister(bam_path = bam_path, which_gr = ir_filt) %>% 
    dplyr::mutate(sample_id = bam_name, 
                  arm = ifelse(strand == "+", "antisense", "sense"))
  
}) %>% 
  dplyr::bind_rows(.)

### plot as separate samples
# split data.frame to list
plot_list <- 
  plot_facet_df %>% 
  split(., plot_facet_df$sample_id)

# plot in a loop
purrr::map(names(plot_list), function(bam_name){
  
  # print bam name
  cat(bam_name, "\n")
  
  # create table for plot
  plot_df <- 
    plot_list[[bam_name]] %>% 
    dplyr::mutate(read_register = replace(read_register, read_register == 0, 22)) %>% 
    dplyr::arrange(read_register)
  
  # duplicate start and end
  plot_df <- 
    dplyr::bind_rows(plot_df %>% dplyr::filter(strand == "+"), 
                     plot_df %>% dplyr::filter(strand == "+", read_register == 1),
                     plot_df %>% dplyr::filter(strand == "-"), 
                     plot_df %>% dplyr::filter(strand == "-", read_register == 1)) %>% 
    dplyr::mutate(read_register = str_c("r", read_register))
  
  # create plot
  p <- 
    plot_ly(plot_df, 
            type = "scatterpolar", 
            mode = "lines",
            r = ~percent, 
            theta = ~read_register, 
            fill = "toself",
            color = ~arm, 
            colors = c("#00BFC4", "#F8766D"), 
            fillcolor = rep(c("#00BFC4", "#F8766D"), each = 23)) %>% 
    layout(polar = list(angularaxis = list(rotation = 90, direction = "clockwise", 
                                           tickmode = "array", 
                                           tickvals = plot_df$read_register, 
                                           ticktext = str_remove(plot_df$read_register, "r"), 
                                           tickfont = list(size = 19), 
                                           showticklabels = F, 
                                           gridcolor = toRGB("gray60"), 
                                           linecolor = toRGB("gray60"), 
                                           ticks = "", 
                                           gridwidth = 2.5,
                                           linewidth = 2.5), 
                        radialaxis = list(visible = T, range = c(0, 100),
                                          tickmode = "array", 
                                          tickvals = c(0, 25, 50, 75, 100),
                                          tickprefix = "     ", 
                                          ticks = "",
                                          tickfont = list(size = 19),
                                          showticklabels = FALSE, 
                                          gridcolor = toRGB("gray60"),
                                          showline = FALSE, 
                                          gridwidth = 2.5)), 
           showlegend = FALSE)
  
  export(p, file.path(outpath, str_c(bam_name, "IR_phase.png", sep = ".")))
  
  ### black and white
  p <- 
    plot_ly(plot_df, 
            type = "scatterpolar", 
            mode = "lines",
            r = ~percent, 
            theta = ~read_register, 
            fill = "toself",
            color = ~arm,
            colors = c(gplots::col2hex("grey50"), gplots::col2hex("black")),
            fillcolor = rep(c(gplots::col2hex("grey50"), gplots::col2hex("black")), each = 23)) %>% 
    layout(polar = list(angularaxis = list(rotation = 90, direction = "clockwise", 
                                           tickmode = "array", 
                                           tickvals = plot_df$read_register, 
                                           ticktext = str_remove(plot_df$read_register, "r"), 
                                           tickfont = list(size = 19), 
                                           showticklabels = F, 
                                           gridcolor = toRGB("gray90"), 
                                           linecolor = toRGB("gray90"), 
                                           ticks = "", 
                                           gridwidth = 2.5,
                                           linewidth = 2.5), 
                        radialaxis = list(visible = T, range = c(0, 100),
                                          tickmode = "array", 
                                          tickvals = c(0, 25, 50, 75, 100),
                                          tickprefix = "     ", 
                                          ticks = "",
                                          tickfont = list(size = 19),
                                          showticklabels = FALSE, 
                                          gridcolor = toRGB("gray90"),
                                          showline = FALSE, 
                                          gridwidth = 2.5)), 
           showlegend = FALSE)
  
  export(p, file.path(outpath, "bw", str_c(bam_name, "IR_phase.bw.png", sep = ".")))
  
  # return
  return(bam_name)
  
})