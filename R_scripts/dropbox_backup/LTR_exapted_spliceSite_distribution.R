# Hadley-verse
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)

# Bioconducor libraries
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(Biostrings)
library(RWebLogo)

options(bitmapType = "cairo")

changeWD <- function(wd){
  if(!dir.exists(wd)){
    dir.create(wd, showWarnings = TRUE, recursive = FALSE, mode = "0755")
    setwd(wd)
    cat("Working directory created and set to", getwd(), "\n")
  } else {
    setwd(wd)
    cat("Working directory set to", getwd(), "\n")
  }
}

changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/polyA_splice_pattern_distribution/exapted_LTRs")

################################################################################## reading data
# Vedran's figure 3 data 
repeats_list <- readRDS("/common/WORK/vfranke/Projects/PSvoboda_MT/Results/MT_FindRepeats/mm/160907.mm.Results.rds", refhook = NULL)

# LTR order 
LTR_order <- c("ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data with consensus length (Maja)
consensus_length <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/LTR_alignment/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB")) %>% 
  group_by(LTR_class) %>%
  summarize(consensus_width = max(consensus_width)) %>% 
  dplyr::filter(!str_detect(LTR_class, "MLT")) %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order)) 

################################################################################## read annotated exons' splice donor sites
# read ENSEMBL and UCSC splice donors
splice_donors <- 
  rbind(read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/splicing/Ensembl_GRCm38.86.20161128_splice_donor.txt", delim =  "\t"), 
        read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/splicing/UCSC_knownGene_mm10_20161126_splice_donor.txt", delim =  "\t")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

################################################################################## read splice junctions coordinates
# get splice junction starts in GV
all_strands <- setNames(c("*", "+", "-"), as.character(0:2))
all_intron_motifs <- setNames(c("non-canonical", "GT/AG", "CT/AC", "GC/AG", "CT/GC", "AT/AC", "GT/AT"), as.character(0:6))

sj_splice_start <- 
  read_delim("/common/WORK/vfranke/Projects/Oocyte_Transcriptome_07102012/Data/InHouseData/Fugaku_RNASeq/Mapped/STAR_PairedEnd_mm10_v2/s_GV.WE/s_GV.WESJ.out.tab", 
             delim = "\t", col_names = c("seqnames", "start", "end", "strand", "intron_motif", "annot", "n_uniq", "n_multi", "overhang")) %>%
  mutate(strand = all_strands[as.character(strand)], 
         intron_motif = all_intron_motifs[as.character(intron_motif)], 
         SJ_fullName = paste0(seqnames, ":",
                              start, "-", 
                              end, ":", 
                              strand)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
  GenomicRanges::resize(., fix = "start", width = 1)

################################################################################## 
################################################################################## 
# full contribution to a 5’ exon in protein-coding genes (red in figure 3D)
selected_LTRs <- 
  repeats_list$RepeatsSelected_MALR %>% 
  filter(ex.category == "5' exon", 
#          ex.category.complete == "Complete",
#          gene_biotype == "protein_coding",
         !is.na(repClassCust)) %>% 
  dplyr::select(coordinates = rep.coord, strand = rep.strand, LTR_subclass = repName, LTR_class = repClassCust) %>% 
  left_join(consensus_length, by = "LTR_class") %>% 
  mutate(fullName = paste0(coordinates, ":", 
                           strand, "|", 
                           LTR_subclass, "|",
                           LTR_class)) %>% 
  separate(coordinates, c("seqnames", "coordinates"), ":") %>% 
  separate(coordinates, c("start", "end"), "-") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# plot name
# plot_name <- "LTR_exapted_fullContribution_proteinCoding"
# plot_name <- "LTR_exapted_fullContribution_all"
plot_name <- "LTR_exapted_fullAndPartialContribution_all"

################################################################################## get splice start distribution
# overlap splice junction starts with selected LTRs
sj_splice_start_LTR <- sj_splice_start[subjectHits(findOverlaps(selected_LTRs, sj_splice_start))] 
sj_splice_start_LTR$fullName <- selected_LTRs[queryHits(findOverlaps(selected_LTRs, sj_splice_start))]$fullName

# overlap splice starts with annotated splice donors, filter out those overlaping 
exon_splice_donors <- sj_splice_start_LTR[queryHits(findOverlaps(sj_splice_start_LTR, splice_donors))]$SJ_fullName
sj_splice_start_LTR <- 
  sj_splice_start_LTR[!(sj_splice_start_LTR$SJ_fullName %in% exon_splice_donors)] %>% 
  as.data.frame() %>% 
  dplyr::select(splice_start = start, splice_motif = intron_motif, fullName)

# join with LTRs
splice_start_distribution <- 
  as.data.frame(selected_LTRs) %>% 
  inner_join(sj_splice_start_LTR) %>% 
  mutate(relative_splice_start = ifelse(strand == "+", splice_start - start, end - splice_start)) %>% 
  dplyr::select(start = relative_splice_start, names = fullName, consensus_width) %>% 
  dplyr::filter(start < consensus_width) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order), 
         pattern = "splice")

################################################################################## get polyA pattern distribution
# get LTRs which have splice donor splice site for polyA pattern matching
selected_LTRs_splice <- selected_LTRs[selected_LTRs$fullName %in% sj_splice_start_LTR$fullName]

# get sequences 
LTR_seqences <- 
  getSeq(x = Mmusculus, selected_LTRs_splice) %>% 
  setNames(selected_LTRs_splice$fullName)

# get pattern start
polyA_distribution <- 
  unlist(vmatchPattern("AATAAA", LTR_seqences, fixed = F)) %>% 
  as.data.frame() %>%
  dplyr::select(start, names) %>% 
  separate(names, c("LTR_range", "LTR_subclass", "LTR_class"), sep = "\\|") %>% 
  left_join(consensus_length, by = "LTR_class") %>% 
  dplyr::filter(start < consensus_width) %>% 
  mutate(LTR_class = factor(LTR_class, levels = LTR_order), 
         pattern = "polyA")

################################################################################## plot
# join polyA/splice distribution
distribution_all <-
  rbind(splice_start_distribution, polyA_distribution) %>% 
  mutate(start = -(consensus_width - start))

consensus_length_negative <- 
  consensus_length %>% 
  mutate(consensus_width = -consensus_width)

# plot pattern start distribution as jitterplot
ggplot(data = distribution_all, aes(x = LTR_class, start)) +
  geom_bar(data = consensus_length_negative, aes(x = LTR_class, y = consensus_width), stat = "identity", alpha = 0.2) +
  geom_jitter(aes(colour = pattern), size = 0.1, width = 0.25) +
  xlab("LTR class") +
  ylab("pattern start position") +
  ylim(c(min(consensus_length_negative$consensus_width), 0)) +
  scale_colour_manual(values = c("black", "red"), 
                      breaks = c("polyA", "splice")) +
  scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggsave(paste0(plot_name, ".png"))
