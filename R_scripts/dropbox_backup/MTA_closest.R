# repeatMasker
rptmsk <- 
  read_delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20161012.txt.gz", delim =  "\t") %>% 
  dplyr::rename(LTR_subclass = element_name) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

selected_LTRs <- 
  repeats_list$RepeatsSelected_MALR %>% 
  filter(ex.category == "5' exon",
         ex.category.complete == "Complete", 
         !is.na(repClassCust)) %>% 
  dplyr::select(coordinates = rep.coord, strand = rep.strand, LTR_subclass = repName, LTR_class = repClassCust) %>% 
  mutate(fullName = paste0(coordinates, ":", 
                           strand, "|", 
                           LTR_subclass, "|",
                           LTR_class)) %>% 
  dplyr::filter(LTR_class == "MTA") %>% 
  separate(coordinates, c("seqnames", "coordinates"), ":") %>% 
  separate(coordinates, c("start", "end"), "-") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# find the closest one
selected_LTRs_closest <- 
  rptmsk[nearest(selected_LTRs, rptmsk, ignore.strand = FALSE)] %>% 
  as.data.frame() %>% 
  dplyr::select(closest = LTR_subclass, closest_start = start, closest_end = end)

# join with LTRs
selected_LTRs_filtered <- 
  selected_LTRs %>% 
  as.data.frame() %>% 
  cbind(selected_LTRs_closest) %>% 
  dplyr::filter(closest == "MTA_Mm-int") %>% 
  write_csv("exapted_MTA_with_MTAMmInt.csv")

