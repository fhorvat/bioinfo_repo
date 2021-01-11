### INFO: 
### DATE: Wed Jul 15 10:34:22 2020
### AUTHOR: @pepap 20200714, modified by Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/substitution_rate/random_sampled_LTRs.all_classes")

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

library(Biostrings)
library(parallel)
library(DECIPHER)
library(ape)
library(phangorn)

######################################################## SOURCE FILES

######################################################## FUNCTIONS
# (1) DNAStringSet => multi-alignment
pep.alnSeqs <- function(inpDSS, ITER = 100, REFIN = 100, PROC = 1, VERBOSE = T){
  
  # Align sequences
  xALN <- DECIPHER::AlignSeqs(myXStringSet = inpDSS, iterations = ITER, refinements = REFIN, processors = PROC, verbose = VERBOSE)
  
  return(xALN)
  
}

# (2) Multi-alignment => substitution matrix
pep.mutRate <- function(inpDALN, INDELS = F, NDIGITS = 5){
  
  # Create output matrix
  out.mat  <- matrix(nrow = length(inpDALN), ncol = length(inpDALN), dimnames = list(names(inpDALN), names(inpDALN)))
  
  # Convert alignments to character strings
  list_CHR <- strsplit(x = as.character(inpDALN), split = "", fixed = T)
  
  # Fill matrix by pair-sequence comparison
  for(i in seq(from = 1, to = length(inpDALN))) {
    
    for(j in seq(from = 1, to = length(inpDALN))) {
      
      # Filter out alignment-positions which equal "-" in both sequences
      not_dash  <- ((list_CHR[[i]] != "-") | (list_CHR[[j]] != "-"))
      
      # Filter out matches
      N_matches <- sum(list_CHR[[i]][not_dash] == list_CHR[[j]][not_dash])
      
      if(INDELS){
        
        # Filter out mismatches + indels
        N_substit <- sum((list_CHR[[i]][not_dash] != list_CHR[[j]][not_dash]))
        
      }else{
        
        # Filter out mismatches & remove indels
        N_substit <- sum((list_CHR[[i]][not_dash] != list_CHR[[j]][not_dash]) & (list_CHR[[i]][not_dash] != "-") & (list_CHR[[j]][not_dash] != "-"))
        
      }
      
      # Calculate substitution rate
      x_val <- round(x = (N_substit / (N_matches + N_substit)), digits = NDIGITS)
      out.mat[i, j] <- x_val
      out.mat[j, i] <- x_val
      
      # if N matches == 0, substitution rate is 1
      if(is.nan(x_val)){
        out.mat[i, j] <- 1
        out.mat[j, i] <- 1
      }
    }
    
  }
  
  # return
  return(out.mat)
  
}

# (3) Substitution matrix => shortest edge-lengths
pep.phylo <- function(inpMUTMAT, RETURN.PHYLO = F){
  
  # Transform the substitution-matrix into phylogenic-tree object
  out.upgma <- phangorn::upgma(D = inpMUTMAT)
  
  if(RETURN.PHYLO){
    return(out.upgma)
  }
  
  # Extract shortest edge-lengths
  out.dist <- out.upgma$edge.length[out.upgma$edge[ , 2] <= nrow(inpMUTMAT)]
  names(out.dist) <- out.upgma$tip.label[out.upgma$edge[out.upgma$edge[, 2] <= nrow(inpMUTMAT), 2]]
  
  return(out.dist)
  
}


# (4) Function performing all previous transformations step-by-step
# n.iter       ... number of iterations within the multi-alignment procedure (default 100)
# n.refin      ... number of refinement steps within the multi-alignment procedure (default 100)
# ncpus        ... number of processors used within the multi-alignment procedure (default 1)
# quiet        ... turn off verbose output from the multi-alignment procedure (default FALSE : if TRUE, the alignment convergence is not reported)
# add.indels   ... include indels into the substitution rate (default FALSE: if TRUE, try to refine your input sequences & remove all additional inserts)
# n.digits     ... number of digits to which the value of substitution rate is rounded (default 5)
# return.phylo ... instead of set of shortest-length values return the original phylogenetic-tree object (default FALSE: if TRUE, you can directly plot the object)

pep.seq2subsRate <- function(inpDNAStringSet, n.iter = 100, n.refin = 100, ncpus = 1, quiet = F, add.indels = F, n.digits = 5, return.phylo = F){
  
  # DNAStringSet => multi-alignment
  cat("\n ++ Input DNAStringSet => multi-alignment ++\n", sep = "")
  out.aln <- pep.alnSeqs(inpDSS = inpDNAStringSet, ITER = n.iter, REFIN = n.refin, PROC = ncpus, VERBOSE = (!quiet))
  
  # Multi-alignment => substitution matrix
  cat(" ++ Multi-alignment => substitution matrix ++\n",sep="")
  out.mut <- pep.mutRate(inpDALN = out.aln, INDELS = add.indels, NDIGITS = n.digits)
  
  # Substitution matrix => shortest edge-lengths
  cat(" ++ Substitution matrix => shortest edge-lengths ++\n",sep="")
  out.dis <- pep.phylo(inpMUTMAT = out.mut, RETURN.PHYLO = return.phylo)
  
  # done, return
  cat(" ++ Finished ++\n\n", sep = "")
  return(out.dis)
  
}

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# table with LTRs info path
ltr_tb_path <- file.path(inpath, "LTRs.random_50_per_repName.raw_rmsk.csv")

# list subfamily fasta files
ltr_seq_path <- list.files(inpath, pattern = ".*\\.single\\.fasta", full.names = T)

######################################################## READ DATA
# read table with LTRs info
ltr_tb <- readr::read_csv(ltr_tb_path)

# read individual LTR subfamily sequences
ltr_seq_list <- 
  purrr::map(ltr_seq_path, Biostrings::readDNAStringSet) %>% 
  set_names(., ltr_seq_path %>% basename(.) %>% str_remove(., "LTRs\\.random_50_per_repName\\.") %>% str_remove(., "\\.raw_rmsk\\.single\\.fasta"))

# remove empty sequences
ltr_seq_list <- ltr_seq_list[purrr::map(ltr_seq_list, length) > 0]

######################################################## MAIN CODE
# get substitution rate in LTRs by running functions above
sub_rate <- purrr::map(ltr_seq_list, pep.seq2subsRate, ncpus = 12)

# set names
names(sub_rate) <- names(ltr_seq_list)


### plot
# prepare table for plot
sub_rate_tb <- purrr::map(names(sub_rate), function(ltr_name){
  
  # create table with substitution rates and ltr subfamily
  tibble(sub_rate = sub_rate[[ltr_name]], 
         repSubfamily = ltr_name)
  
}) %>% 
  bind_rows(.) %>% 
  dplyr::arrange(repSubfamily)

# data statistics
sub_rate_stats <- 
  sub_rate_tb %>% 
  group_by(repSubfamily) %>% 
  summarize(median = median(sub_rate), 
            sd = sd(sub_rate), 
            length = length(sub_rate)) %>% 
  dplyr::mutate(SEM = sd / sqrt(length))

# labels for plot
labels_tb <- 
  sub_rate_tb %>% 
  dplyr::group_by(repSubfamily) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::filter(count > 0)
  
# plot
sub_rate_boxplot <- 
  ggplot() +
  stat_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), geom = "errorbar", ) +
  geom_boxplot(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), outlier.colour = NULL, outlier.shape = NA) +
  geom_jitter(data = sub_rate_tb, aes(x = repSubfamily, y = sub_rate), alpha = 0.4, 
              colour = "black", size = 1, width = 0.1, height = 0, show.legend = F) +
  scale_x_discrete(labels = str_c(labels_tb$repSubfamily, 
                                  " (",
                                  labels_tb$count, 
                                  ")"), 
                   drop = TRUE) + 
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())

# save plot
ggsave(plot = sub_rate_boxplot, 
       filename = file.path(outpath, str_c("LTRs.random_50_per_repName", "substitution_rate", "raw_rmsk", "boxplot.png", sep = ".")), 
       width = 30, 
       height = 10)

