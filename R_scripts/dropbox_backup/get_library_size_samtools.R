## get library size in million of reads (number of uniquely mapped reads)
# single-end: samtools view -F 0x904 -c *.bam
# paired-end: samtools view -F 0x904 *.bam | cut -f 1 | sort | uniq | wc -l

library_size <- sapply(X = filenames, FUN = function(X){
  system(command = paste0("samtools view -F 0x904 -c ", X), intern = T)
})

library_size_df <- 
  data.frame(sample = names(library_size), 
             library_size = as.integer(library_size)) %>% 
  mutate(sample = str_replace_all(sample, "\\/.*\\/|.bam", ""), 
         library_size = library_size / 10^6) %>% 
  dplyr::slice(match(c("s_GV", "s_MII", "s_4c", "s_8c", "s_16c", "s_Blast"), sample)) %T>% 
  write_delim(path = "Graf_2014_library_size.txt", delim = "\t", col_names = T)