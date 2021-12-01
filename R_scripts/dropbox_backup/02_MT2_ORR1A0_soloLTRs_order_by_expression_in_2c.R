library(GenomicRanges)
library(CoverageView)
library(Biostrings)
library(seqinr)

library(rtracklayer)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)

library(readr)
library(purrr)
library(dplyr) 
library(stringr)
library(magrittr)
library(tidyr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/Park_2013_mm9")

################################################################################## reading data
# lift over chain to mm9
mm9_chain <- import.chain("/common/WORK/fhorvat/reference/mouse/mm9/mm10ToMm9.over.chain")

# making TxDb object from knownGene gtf from UCSC
knownGenes_gtf_gr <- 
  makeTxDbFromGFF("/common/WORK/fhorvat/reference/mouse/mm10/UCSC_knownGene_mm10_20161126.gtf.gz") %>% 
  genes(.)

# combining both MT2 and ORR1A0 solo LTRs, making GRanges
MT2_ORR1A0_solo_gr <- 
  rbind(read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/MT2_solo_LTRs_Maja_20160803.txt", delim = "\t"), 
                         read_delim("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/ORR1A0_solo_LTRs_Maja_20160803.txt", delim = "\t")) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# finding overlaps between MT2/ORR1A0 solo LTRs and knownGene table, filtering those which overlap
MT2_ORR1A0_solo_knownGenes_overlaps <- findOverlaps(MT2_ORR1A0_solo_gr, knownGenes_gtf_gr)
MT2_ORR1A0_solo_filtered <- 
  MT2_ORR1A0_solo_gr[-queryHits(MT2_ORR1A0_solo_knownGenes_overlaps), ] %>% 
  liftOver(., mm9_chain) %>% 
  unlist()

################################################################################## counting overlaps
# .bam files path
filenames <- file.path(c("/common/DB/vfranke/Base/AccessoryData/Park_2013_GenDev/Merged/s_Oo.SE/s_Oo.SE.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Park_2013_GenDev/Merged/s_2C.SE/s_2C.SE.bam"))

# # library size
# samtools view -F 0x904 -c s_Oo.SE.bam
number_of_reads <- 
  read_lines("Park_2013_library_size.txt") %>% 
  as.integer() %>% 
  divide_by(10^6) %>% 
  set_names(c("s_Oo", "s_2C"))

# counts 
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
se <- summarizeOverlaps(features = MT2_ORR1A0_solo_filtered, 
                        reads = bamfiles["s_2C.SE.bam"], 
                        mode = "Union", 
                        singleEnd = TRUE, 
                        ignore.strand = TRUE)

# calculating FPKM from counts, order by FPKM, get unique ordered ID
MT2_ORR1A0_solo_fpkm_ordered <- 
  as.data.frame(assay(se)) %>% 
  set_colnames("s_2C_fpkm")  %>% 
  mutate(width = width(MT2_ORR1A0_solo_filtered), 
         s_2C_fpkm = s_2C_fpkm / (number_of_reads["s_2C"] * (width / 1000))) %>%
  dplyr::select(-width) %>% 
  cbind(as.data.frame(MT2_ORR1A0_solo_filtered), .) %>% 
  dplyr::arrange(desc(s_2C_fpkm)) %T>% 
  write_delim(path = "MT2_ORR1A0_solo_LTRs_orderedByExpressionIn2cell_Park_mm9.txt", delim = "\t", col_names = T)
