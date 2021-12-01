### INFO: 
### DATE: Tue Sep 11 09:14:41 2018
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
setwd("/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/phasing/alignment_plots")

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
library(ggbio)

######################################################## SOURCE FILES
lib_path <- "/common/WORK/fhorvat/code_library/R_scripts"
source(file.path(lib_path, "wideScreen.R"))
source(file.path(lib_path, "headt.R"))
source(file.path(lib_path, "asdf.R"))

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# get .BED with inverted repeat coordinates
ir_path <- "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/Analysis/Documentation/plasmid_sequences/all_plasmids.bed"

######################################################## READ DATA
# read coordinates of inverted repeats
ir_gr <- rtracklayer::import.bed(con = ir_path)
start(mosir_gr) <- start(mosir_gr) - 1

######################################################## MAIN CODE
# experiment paths
experiment_path <- 
  "/common/WORK/fhorvat/Projekti/Svoboda/smallRNASeq_2018/NIH3T3_transfected.2018/Data/Mapped/STAR_mm10.all_plasmids" %>% 
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
for(bam_path in bam_paths){
  
  # get bam name
  bam_name <- basename(bam_path) %>% str_remove(., ".bam")
  
  # get IR name
  ir_name <- str_remove_all(bam_name, "^s_|_r[1,2].SE.mis_[0-2]") %>% str_replace(., "pCag_EGFP", "pCag-EGFP")
  
  # subset IR
  ir_subset <- ir_gr[seqnames(ir_gr) == ir_name]
  
  # print bam name
  cat(bam_name, "\n")
  
  # read bam
  bam_gr <- 
    GenomicAlignments::readGAlignmentsList(file = bam_path, 
                                           use.names = TRUE, 
                                           param = ScanBamParam(which = ir_subset, 
                                                                flag = scanBamFlag(isMinusStrand = F))) %>% 
    unlist(.)
  
  # take only reads between 21-23nt
  bam_gr <- bam_gr[str_detect(cigar(bam_gr), "21M|22M|23M")]
  
  # add overlap with arm
  values(bam_gr)$arm <- subjectHits(findOverlaps(bam_gr, ir_subset))
  
  # normalize start of IR
  bam_norm <- 
    bam_gr %>% 
    GRanges(.) %>% 
    as.data.frame(.) %>%
    tibble::as.tibble(.) %>%
    dplyr::select(seqnames, start, end, arm) %>% 
    dplyr::mutate(start_norm = ifelse(test = (arm == 1), 
                                      yes = (start - start(ir_subset[1]) + 1), 
                                      no = (end(ir_subset[2]) - end + 1)), 
                  end_norm = ifelse(test = (arm == 1), 
                                    yes = (end - start(ir_subset[1]) + 1),
                                    no = (end(ir_subset[2]) - start + 1))) %>% 
    dplyr::select(seqnames, start = start_norm, end = end_norm, arm) %>% 
    dplyr::filter(start > 0)
  
  ### collapsed reads
  # get GRanges
  bam_collapsed <- 
    bam_norm %>%
    dplyr::group_by(seqnames, start, end, arm) %>%
    dplyr::summarise(count = n()) %>%
    GRanges(.)
  
  # arm 1
  reads_plot_arm1 <- 
    ggbio::autoplot(bam_collapsed[bam_collapsed$arm == 1], aes(fill = count, color = count), size = 0.35) +
    scale_fill_continuous(type = "viridis", limits = c(1, max(bam_collapsed$count))) +
    scale_color_viridis_c() +
    guides(color = F)
  
  # arm 2
  reads_plot_arm2 <- 
    ggbio::autoplot(bam_collapsed[bam_collapsed$arm == 2], aes(fill = count, color = count), size = 0.35) +
    scale_fill_continuous(type = "viridis", limits = c(0, max(bam_collapsed$count))) +
    scale_color_viridis_c() +
    scale_y_reverse(breaks = NULL, minor_breaks = NULL, labels = NULL) +
    guides(fill = F, color = F)
  
  # create and save ggbio plot
  ggbio::tracks(`arm 1` = reads_plot_arm1,
                `arm 2` = reads_plot_arm2,
                heights = c(1, 1), 
                theme = theme_clear(), 
                label.text.cex = 0.75) +
    ggsave(filename = file.path(outpath, experiment, str_c(bam_name, "collapsed", "IR_alignments.png", sep = ".")), width = 20, height = 10)
  
}




