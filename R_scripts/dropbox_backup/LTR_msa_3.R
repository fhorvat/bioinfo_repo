################################################################################### calculate mean FPKM expression in 3' half of LTR in GV, order 
# all chosen LTRs in repeat masker
rptmsk_LTR <-
  rptmsk %>%
  filter(grepl("LTR", element_class),
         grepl("MLT|MT|ORR", element_name),
         !grepl("int", element_name),
         !grepl("MLTR", element_name)) %>%
  mutate(fullName = paste0(seqnames, ":",
                           start, "-",
                           end, ":",
                           strand, "|",
                           element_name)) %>%
  dplyr::select(-element_class) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

# resize to 3' half
rptmsk_LTR_3prime_half <- GenomicRanges::resize(rptmsk_LTR, width = round(width(rptmsk_LTR) / 2), fix = "end")
rptmsk_LTR_3prime_half_expressed <- rptmsk_LTR_3prime_half[rptmsk_LTR_3prime_half$fullName %in% rptmsk_LTR_expressed$fullName]

register(MulticoreParam())
se_full <- summarizeOverlaps(features = rptmsk_LTR_3prime_half, 
                        reads = BamFileList(sample_df_LTR$track_path, yieldSize = 2000000), 
                        mode = "Union", 
                        singleEnd = F, 
                        ignore.strand = TRUE)

register(MulticoreParam())
se_sub <- summarizeOverlaps(features = rptmsk_LTR_3prime_half_expressed, 
                            reads = BamFileList(sample_df_LTR$track_path, yieldSize = 2000000), 
                            mode = "Union", 
                            singleEnd = F, 
                            ignore.strand = TRUE)

se_sub_df <- 
  as.data.frame(assay(se_sub)) %>% 
  mutate(fullName = rptmsk_LTR_3prime_half_expressed$fullName) 

se_full_df <- 
  as.data.frame(assay(se_full)) %>% 
  mutate(fullName = rptmsk_LTR_3prime_half$fullName)

se_both <- 
  left_join(se_sub_df, se_full_df, by = "fullName") %>% 
  dplyr::select(fullName, subGV = s_GV.WE.bam.x, fullGV = s_GV.WE.bam.y) %>% 
  mutate(sub_full = subGV - fullGV)