# find o
bam_gr_overlaps <- 
  GenomicRanges::findOverlaps(bam_gr, minoverlap = 20) %>% 
  as.data.frame(.) %>% 
  dplyr::filter(queryHits != subjectHits) %>% 
  apply(., 1, sort) %>% 
  t(.) %>% 
  unique(.) %>% 
  as.data.frame(.) %>% 
  dplyr::rename(query = V1, subject = V2) %>% 
  tibble::as_tibble(.)


x <- punion(bam_gr[bam_gr_overlaps$query], bam_gr[bam_gr_overlaps$subject])