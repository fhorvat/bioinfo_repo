### INFO: 
### DATE: Sat Sep 28 23:34:42 2019
### AUTHOR: Filip Horvat
rm(list = ls()); gc()
options(bitmapType = "cairo")
wideScreen()

######################################################## WORKING DIRECTORY
# set working dir
setwd("/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/FLI_elements/IAP/annotation")

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

library(BSgenome.Maur.UCSC.Siomi)
library(seqinr)
library(Biostrings)
library(systemPipeR)
library(stringdist)
library(msa)
library(DECIPHER)
library(seqinr)
library(ape)

######################################################## SOURCE FILES

######################################################## FUNCTIONS

######################################################## PATH VARIABLES
# set inpath 
inpath <- getwd()

# set outpath
outpath <- getwd()

# genome path
genome_dir <- "/common/DB/genome_reference/golden_hamster/Siomi_assembly.fixed"

# insertions table path
insertions_path <- "/common/WORK/fhorvat/Projekti/Svoboda/piRNA.Zuzka/Analysis/2020_paper/retrotransposons_annotation/golden_hamster.Siomi_assembly/LTRs/potentially_young_LTRs"
insertions_path <- file.path(insertions_path, "potentially_young_LTRs.edit_distance.10.ORFs.csv")

######################################################## READ DATA
# read insertion table
insertions_tb <- readr::read_csv(insertions_path)

######################################################## MAIN CODE
# tidy table, filter only IAPs, order by ORF length (GAG = 576 aa, POL = 905 aa)
iap_tb <- 
  insertions_tb %>% 
  dplyr::filter(ltr_class == "IAP") %>% 
  dplyr::select(insertion_coordinates, strand, ltr_name, rmsk_id, longest_orf_1, longest_orf_2, insertion_class) %>% 
  tidyr::separate(insertion_coordinates, into = c("seqnames", "start", "end"), sep = " ") %>% 
  dplyr::mutate(start = as.numeric(start), 
                end = as.numeric(end), 
                width = end - start + 1) %>% 
  dplyr::arrange(desc(longest_orf_1 + longest_orf_2))

# # save as table
# readr::write_csv(iap_tb, file.path(outpath, "IAP.potentially_young.ordered_by_ORFs.20201031.new.csv"))


### top N numbers, based on Petr's assesment 
n_top <- 110

# filter by ORF length 
iap_tb %<>% 
  dplyr::top_n(n_top, longest_orf_1 + longest_orf_2)

# # save as table
# readr::write_csv(iap_tb, file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.top_", n_top, ".csv")))

# to GRanges
iap_gr <-
  iap_tb %>%
  GRanges(.)
names(iap_gr) <- str_c(mcols(iap_gr)$ltr_name, mcols(iap_gr)$rmsk_id, sep = ".")
names(iap_gr) <- str_replace_all(names(iap_gr), "-", "_")

# # save as .bed
# rtracklayer::export.bed(iap_gr, file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.top_", n_top, ".bed")))

# get sequences
iap_seq <- getSeq(BSgenome.Maur.UCSC.Siomi, iap_gr)

# # # save as fasta
# Biostrings::writeXStringSet(x = iap_seq, filepath = file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.top_", n_top, ".fasta")))


### get consensus by repeat names
# split by repeat names
iap_seq_list <- split(iap_seq, iap_gr$ltr_name)

## MSA and get consesus
iap_seq_consensus <- purrr::map(names(iap_seq_list), function(seq_name){
  
  seq_name <- names(iap_seq_list)[1]
  
  # do the MSA
  seq_msa <- msa::msa(iap_seq_list[[seq_name]], method = "ClustalOmega")
  # seq_msa <- DECIPHER::AlignSeqs(myXStringSet = iap_seq_list[[seq_name]], iterations = 100, refinements = 100, processors = 1, verbose = T)
  
  # write as fasta
  seq_msa@unmasked %>%
    Biostrings::writeXStringSet(., file = file.path(outpath, 
                                                    str_c("IAP.potentially_young.ordered_by_ORFs.20201031.", 
                                                          seq_name, 
                                                          ".msa.fasta")))
  
  # get consensus scores
  data(BLOSUM62)
  consensus_score <- msaConservationScore(x = seq_msa, substitutionMatrix = BLOSUM62, gapVsGap = 0,
                                          type = "upperlower", thresh = c(80, 0), ignoreGaps = T)

  ### plot the consensus score
  # create table
  plot_tb <-
    tibble(letter = names(consensus_score),
           score = unname(consensus_score)) %>%
    dplyr::mutate(letter = toupper(letter),
                  position = 1:n())

  # plot as barplot
  barplot_plot <-
    ggplot() +
    geom_density(data = plot_tb[1:1000, ],
                 mapping = aes(x = position, y = score), stat = "identity") +
    # scale_fill_manual(values = c(miRNA = "#70ad47", mRNA = "#ffc000", repeats = "#ff0000", other_mapped = "#000000")) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

  # save
  ggsave(plot = barplot_plot,
         filename = file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.", seq_name, ".top_", n_top, ".msa.score.png")),
         width = 30, height = 5,
         limitsize = F)
  
  # get the consensus
  seq_consensus <-
    msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
    toupper(.) %>%
    Biostrings::DNAStringSet(.)
  
  # write consensus 
  Biostrings::writeXStringSet(seq_consensus, 
                              file = file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.", 
                                                              seq_name, 
                                                              ".consensus.fasta")))
  
  # return
  return(seq_consensus)
  
})

# # join together
# iap_seq_consensus <- do.call(c, iap_seq_consensus)
# names(iap_seq_consensus) <- names(iap_seq_list)
# 
# # save the consensus
# Biostrings::writeXStringSet(iap_seq_consensus, file = file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.top_", n_top, ".consensus.fasta")))


## MSA together
# MSA 
seq_msa <- msa::msa(iap_seq, method = "ClustalOmega")

# write as fasta
seq_msa@unmasked %>%
  Biostrings::writeXStringSet(., file = file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031", "together", "msa.fasta", sep = ".")))

# get the consensus
seq_consensus <-
  msaConsensusSequence(seq_msa, type = "upperlower", thresh = c(80, 0), ignoreGaps = T) %>%
  toupper(.) %>%
  Biostrings::DNAStringSet(.)

# write consensus 
Biostrings::writeXStringSet(seq_consensus, 
                            file = file.path(outpath, str_c("IAP.potentially_young.ordered_by_ORFs.20201031.", 
                                                            "together", 
                                                            ".consensus.fasta")))


# # convert to seqinr format
# msa_seqinr <- msaConvert(seq_msa, type = "seqinr::alignment")
# 
# # calculate distance matrix
# msa_dist <- seqinr::dist.alignment(msa_seqinr, "identity")
# 
# # get tree and plot
# msa_tree <- ape::njs(msa_dist)
# msa_tree <- ape::ladderize(msa_tree)
# 
# # save with labels
# png(filename = file.path(outpath, str_c("IAP.potentially_young.filtered_by_ORFs", "phylo_tree", "DECIPHER.msa", "labels", "png", sep = ".")), width = 2000, height = 2000, units = "px")
# plot.phylo(msa_tree, edge.width = 2, font = 1, show.tip.label = TRUE, cex = 1)
# edgelabels(round(msa_tree$edge.length, 3), bg = "black", col = "white", font = 2)
# axisPhylo(backward = F)
# dev.off()


# ### expand ranges
# # to GRanges
# seqlevels(line1_gr) <- seqlevels(BSgenome.Cgriseus.Ensembl.CriGriPICR)
# seqlengths(line1_gr) <- seqlengths(BSgenome.Cgriseus.Ensembl.CriGriPICR)
# line1_gr_expand <- line1_gr + 10000
# line1_gr_expand <- trim(line1_gr_expand)
# 
# # get sequences
# line1_seq_expand <- getSeq(BSgenome.Cgriseus.Ensembl.CriGriPICR, line1_gr_expand)
# 
# # save as fasta
# Biostrings::writeXStringSet(x = line1_seq_expand, filepath = file.path(outpath, "LINE1.CriGriPICR.5K_to_7k.ORFs.expand_10k.fasta"))
# 
# 
# # save
# saveRDS(line1_gr, file = file.path(outpath, str_c("LINE1_annotation.rmsk.GRanges.RDS")))
# 
# # save as .bed
# rtracklayer::export.bed(line1_gr, file.path(outpath, "LINE1_annotation.rmsk.bed"))
# 
# # save as SAF
# line1_saf <-
#   line1_gr %>%
#   as_tibble(.) %>%
#   dplyr::select(GeneID = rmsk_id,	Chr	= seqnames, Start	= start, End = end, Strand = strand) %>%
#   readr::write_delim(., file.path(outpath, "LINE1_annotation.rmsk.saf"), delim = "\t")

