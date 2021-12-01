# Hadley-verse
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(reshape2)
library(magrittr)
library(splitstackshape) 

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

changeWD("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/pattern_distribution")

################################################################################## reading data
# repeatMasker
rptmsk <- 
  read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t") %>% 
  dplyr::rename(LTR_subclass = element_name)

# LTR order 
LTR_order <- c("MLT1", "MLT2", 
               "ORR1G", "ORR1F", "ORR1E", "ORR1D", "ORR1C", "ORR1B", "ORR1A", 
               "MTE", "MTE2", "MTD", "MTC", "MTB", "MTA", "MT2A", "MT2B", "MT2C", "MT2")

# LTR data
LTR_data <- 
  read_csv("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/LTR_alignment/files/seq_logo/LTR_ERVL_data_table.csv") %>% 
  dplyr::select(LTR_class, LTR_subclass = RepeatMasker_ID, consensus_width = consensus_sequence_length) %>% 
  mutate(LTR_class = replace(LTR_class, LTR_class == "MTB2", "MTB"), 
         LTR_class = factor(LTR_class, levels = LTR_order)) 

# repeatMasker filtered
rptmsk_filter <- 
  right_join(rptmsk, LTR_data, by = "LTR_subclass") %>% 
  mutate(element_width = (end - start + 1)) %>% 
  filter(element_width > (consensus_width - (0.05 * consensus_width)), 
         element_width < (consensus_width + (0.05 * consensus_width))) %>% 
  dplyr::select(seqnames, start, end, strand, LTR_subclass, LTR_class) %>% 
  mutate(fullName = paste0(seqnames, ":",
                           start, "-", 
                           end, ":", 
                           strand, "|", 
                           LTR_subclass, "|", 
                           LTR_class))

################################################################################## 
# get MTA
LTR_random <- 
  rptmsk_filter %>% 
  dplyr::filter(LTR_class %in% c("MTA", "MTB"))

LTR_random <- makeGRangesFromDataFrame(LTR_random, keep.extra.columns = T)
# write_csv(x = as.data.frame(LTR_random), path = "LTR_random200_perClass.csv")

# get sequences 
LTR_seqences <- getSeq(x = Mmusculus, LTR_random)
names(LTR_seqences) <- LTR_random$fullName

################################################################################## read annotated exons' splice donor sites
# read ENSEMBL and UCSC splice donors
splice_donors <- 
  rbind(read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing/Ensembl_GRCm38.86.20161128_splice_donor.txt", delim =  "\t"), 
        read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/splicing/UCSC_knownGene_mm10_20161126_splice_donor.txt", delim =  "\t")) %>%
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
selected_LTRs <- LTR_random

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
  dplyr::select(splice_start, fullName) %>% 
  group_by(fullName) %>% 
  summarise(splice_start = paste(splice_start, collapse = ",")) %>% 
  cSplit(., "splice_start", ",") %>% 
  as.data.frame()  

splice_start_distribution$n_unique <- apply(splice_start_distribution[, 2:8], 1, function(i) length(unique(i[!is.na(i)])))
splice_start_distribution <- splice_start_distribution[splice_start_distribution$n_unique > 1, ]
write_csv(x = splice_start_distribution, path = "LTR_MTA_MTB_withMultiple_spliceSites.csv")
