library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/class_consensus")

################################################################################## functions
classConsensus <- function(LTR_name){
  
  # get all consensus sequences of one class of LTR
  LTR_consensus <- ltr_consensus_all_split[[LTR_name]]
  
  if(length(LTR_consensus) > 1){
    # write sequences as FASTA
    write.fasta(as.list(LTR_consensus),
                nbchar = 80,
                names = names(LTR_consensus),
                as.string = TRUE,
                file.out = paste0(LTR_name, "_consensus.fasta"),
                open = "w")
    
    # get sequences aligned by Clustal Omega
    LTR_seq_aligned <- ClustalO(infile = paste0(getwd(), "/", LTR_name, "_consensus.fasta"), 
                                outfile = paste0(getwd(), "/", LTR_name, "_consensus_aligned.fasta"),
                                threads = 8, 
                                force = T)
    
    # melt data.frame for plot
    LTR_seq_aligned_df <- 
      as.data.frame(as.matrix(LTR_seq_aligned), stringsAsFactors = F) %>% 
      mutate(ID = rownames(.), 
             ID = factor(ID, levels = rev(ID))) %>% 
      melt(id.vars = "ID") %>% 
      mutate(value = replace(value, value == "-", NA), 
             value = factor(value, levels = c("A", "C", "G", "T")))
    
    # plot
    ggplot(LTR_seq_aligned_df, aes(x = variable, y = ID)) +
      geom_tile(aes(fill = value)) +
      scale_fill_manual(values = color_pallete, 
                        breaks = c("A", "C", "G", "T")) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 3, hjust = 1, vjust = 0),
            #           axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      ggsave(paste0(LTR_name, "_consensus_msa.pdf"))
    
    return(paste0("class ", LTR_name, " aligned"))
    
  }else{
    return(paste0("class ", LTR_name, " has only one subclass"))
  }
}

################################################################################## reading data
comboClasses_list <- list(
  "MLT1" = c("MLT1A", "MLT1A0", "MLT1A1", "MLT1B", "MLT1C", "MLT1D", "MLT1E", "MLT1E1", "MLT1E1A", "MLT1E2", "MLT1E3", "MLT1F", 
             "MLT1F1", "MLT1F2", "MLT1G", "MLT1G1", "MLT1G3", "MLT1H", "MLT1H1", "MLT1H2", "MLT1I", "MLT1J", "MLT1J1", "MLT1J2", 
             "MLT1K", "MLT1L", "MLT1M", "MLT1N2", "MLT1O"),
  "MLT2" = c("MLT2B1", "MLT2B2", "MLT2B3", "MLT2B4", "MLT2B5", "MLT2C1", "MLT2C2", "MLT2D", "MLT2E", "MLT2F"),
  "MT2"   = c("MT2_Mm"),
  "MT2A"  = c("MT2A"),
  "MT2B"  = c("MT2B", "MT2B1", "MT2B2"),
  "MT2C"  = c("MT2C_Mm"),
  "MTA"   = c("MTA_Mm"),
  "MTB"   = c("MTB", "MTB_Mm"),
  #  "MTB2"  = c("MTB_Mm"),
  "MTC"   = c("MTC"),
  "MTD"   = c("MTD"),
  "MTE"   = c("MTEa", "MTEb"),
  "MTE2"  = c("MTE2a", "MTE2b"),
  "ORR1A" = c("ORR1A0", "ORR1A1", "ORR1A2", "ORR1A3", "ORR1A4"),
  "ORR1B" = c("ORR1B1", "ORR1B2"),
  "ORR1C" = c("ORR1C1", "ORR1C2"),
  "ORR1D" = c("ORR1D1", "ORR1D2"),
  "ORR1E" = c("ORR1E"),
  "ORR1F" = c("ORR1F"),
  "ORR1G" = c("ORR1G")
)
comboClasses <- unname(unlist(comboClasses_list))
names(comboClasses) <- comboClasses

# get consensus sequences of all LTRs from repbase
ltr_consensus_all <- lapply(comboClasses, getLTRConsensus)
ltr_consensus_all_split <- split(ltr_consensus_all, rep(names(comboClasses_list), times = unname(sapply(comboClasses_list, length))))

# align and plot
lapply(names(ltr_consensus_all_split), classConsensus)


