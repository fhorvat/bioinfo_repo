library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(CoverageView)
library(Biostrings)
library(seqinr)
library(rtracklayer)

library(readr)
library(purrr)
library(dplyr) 
library(stringr)
library(magrittr)
library(tidyr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/review/coverage_scanning/Park_2013_mm9")

################################################################################## reading data
# loading object
load("/common/WORK/fhorvat/Projekti/Svoboda/Analyses/GenomeResearch_2016/scanning_coverage/gr_mis.Robj")

# lift over chain to mm9
mm9_chain <- import.chain("/common/WORK/fhorvat/reference/mouse/mm9/mm10ToMm9.over.chain")

# filtering for MT2 elements
rmsk_Maja <- 
  gr_mis %>% 
  as.data.frame(.) %>% 
  dplyr::select(seqnames, start, end, strand, element_name = repName, id) %>% 
  dplyr::filter((str_detect(element_name, "MT2")) & (id != 0))

# getting full MT2 elements IDs
fullMT2_id <-  
  rmsk_Maja %>% 
  dplyr::select(start, end, id) %>% 
  group_by(id) %>% 
  nest() %>% 
  mutate(element_width = data %>% 
           map(function(x){
             element_range <- range(x$start, x$end)
             element_range <- element_range[2] - element_range[1]
           }))  %>% 
  dplyr::select(id, element_width) %>% 
  tidyr::unnest() %>% 
  dplyr::filter((element_width > 5900) & (element_width < 6600)) %$% 
  id

# getting full MT2 elements coordinates
MT2full_coord <- 
  rmsk_Maja %>% 
  dplyr::filter(is_in(id, fullMT2_id)) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>% 
  sortSeqlevels(.) %>% 
  sort(.) %>% 
  as.data.frame(.) %>% 
  mutate(end_lead = lead(end), 
         id_lead = lead(id), 
         distance = end_lead - start) %>% 
  dplyr::filter((distance < 7000) & (distance > 0)) %>% 
  dplyr::filter(id == id_lead) %>% 
  dplyr::select(seqnames, start, end = end_lead, strand, element_name, id) %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>% 
  liftOver(., mm9_chain) %>% 
  unlist() %>% 
  magrittr::extract(., width(.) > 6000)
  
# getting coordinates of LTRs from full length MT2 elements
MT2ltr_coord <-
  rmsk_Maja %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>% 
  magrittr::extract(., is_in(.$id, MT2full_coord$id)) %>% 
  liftOver(., mm9_chain) %>% 
  unlist()

################################################################################## counting overlaps
# .bam files path
filenames <- file.path(c("/common/DB/vfranke/Base/AccessoryData/Park_2013_GenDev/Merged/s_Oo.SE/s_Oo.SE.bam",
                         "/common/DB/vfranke/Base/AccessoryData/Park_2013_GenDev/Merged/s_2C.SE/s_2C.SE.bam"))

# library size
number_of_reads <- 
  read_lines("Park_2013_library_size.txt") %>% 
  as.integer() %>% 
  divide_by(10^6) %>% 
  set_names(c("s_Oo", "s_2C"))
  
# counts 
bamfiles <- BamFileList(filenames, yieldSize = 2000000)
se <- summarizeOverlaps(features = MT2ltr_coord, 
                        reads = bamfiles["s_2C.SE.bam"], 
                        mode = "Union", 
                        singleEnd = TRUE, 
                        ignore.strand = TRUE)

################################################################################## 
# calculating FPKM from counts, order by FPKM
MT2ltr_fpkm_ordered_id <- 
  as.data.frame(assay(se)) %>% 
  set_colnames("s_2C_fpkm")  %>% 
  mutate(width = width(MT2ltr_coord), 
         s_2C_fpkm = s_2C_fpkm / (number_of_reads["s_2C"] * (width / 1000)), 
         id = mcols(MT2ltr_coord)$"id") %>% 
  dplyr::select(-width) %>% 
  dplyr::arrange(desc(s_2C_fpkm)) %>% 
  dplyr::filter(!duplicated(id))

# getting ranges of full MT2 ordered by expression of LTRs
MT2full_ordered <- 
  as.data.frame(MT2full_coord) %>%
  left_join(., MT2ltr_fpkm_ordered_id, by = "id") %>% 
  dplyr::arrange(desc(s_2C_fpkm)) %T>% 
  write_delim(path = "MT2_full_orderedByExpressionIn2cell_Park_mm9.txt", delim = "\t", col_names = T)



