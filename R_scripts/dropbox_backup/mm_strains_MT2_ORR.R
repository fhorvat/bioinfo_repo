library("dplyr")
library("readr")
library("ggplot2")
library("GenomicRanges")
library("BSgenome.Mmusculus.UCSC.mm10")
library("seqinr")

# reading mouse repeatMasker table
mm10_rmsk <- read.delim("/common/WORK/fhorvat/reference/mm10/UCSC_repeatMasker_mm10_20160823.fa.out.txt.gz", stringsAsFactors = F, header = T)

# filtering LTRs
mm10_LTRs <- 
  mm10_rmsk %>%
  filter(grepl("LTR", element_class)) %>%
#   mutate(width = end - start + 1) %>%
  mutate(strand = replace(strand, strand == "C", "*")) %>%
  select(-VIZ_ID)

# reading list of MT2 with ORR insertions
MT2_ORR <- 
  read.csv("/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/MT2_ORR_insertion/MT2_fullLengthGVExpressedWithORR1A3int.csv", stringsAsFactors = F, row.names = 1) %>%
  select(-width, -id) %>%
  rename(repName = "element_name")

# creating GRanges and finding overlaps
mm10_LTRs_gr <- makeGRangesFromDataFrame(mm10_LTRs, keep.extra.columns = T)
MT2_ORR_gr <- makeGRangesFromDataFrame(MT2_ORR, keep.extra.columns = T)

# findOverlaps
MT2_ORR_overlaps <- findOverlaps(MT2_ORR_gr, mm10_LTRs_gr)
mm10_MT_ORR <- mm10_LTRs_gr[subjectHits(MT2_ORR_overlaps)]
mm10_MT_ORR <- as.data.frame(mm10_MT_ORR[grepl("ORR", mm10_MT_ORR$element_name)])

# write fasta
ranges_MT_seq <- getSeq(x = Mmusculus,
                        names = mm10_MT_ORR$seqnames,
                        start = as.integer(mm10_MT_ORR$start),
                        end = as.integer(mm10_MT_ORR$end))
ranges_MT_seq <- as.character(ranges_MT_seq)
write.fasta(as.list(ranges_MT_seq),
            nbchar = 80,
            names = paste0(mm10_MT_ORR$seqnames, ":",
                           mm10_MT_ORR$start, "-",
                           mm10_MT_ORR$end, ":",
                           mm10_MT_ORR$strand, "|",
                           mm10_MT_ORR$element_name),
            as.string = TRUE,
            file.out = "/common/WORK/fhorvat/Projekti/Svoboda/Prague_November2016/MT2_ORR_insertion/ORR1A3_int_in_fullLength_MT2.fasta",
            open = "w")


# filtering MT2 LTRs and MERVL-int
mm10_ORR1 <- 
  mm10_LTRs %>%
  filter(grepl("ORR1A3-int", element_name)) 

ggplot(data = mm10_ORR1, aes(width)) +
  geom_histogram(binwidth = 10, aes(fill = element_name)) + 
#   scale_x_continuous(limits = c(300, 400)) +
  facet_grid(element_name ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Mouse strains
mm_path <- list.files("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/annotation/mm_strains/", full.names = T, pattern = "*.csv")
mm_list <- lapply(X = mm_path, FUN = function(x) read.csv(file = x, header = T, stringsAsFactors = F))
names(mm_list) <- gsub("\\/.*\\/|\\..*|_rmsk_LTR", "", mm_path)

mm_MT2_MERVL_list <- lapply(X = mm_list, FUN = function(df){
  df <- df %>% 
    filter(grepl("MT2|MERVL", element_name)) %>%
    mutate(width = end - start + 1) %>%
    mutate(strand = replace(strand, strand == "C", "*")) 
  return(df)
})

mm_MT2_MERVL_df <- do.call(rbind, mm_MT2_MERVL_list)
mm_MT2_MERVL_df$strain <- gsub("\\..*", "", rownames(mm_MT2_MERVL_df))
rownames(mm_MT2_MERVL_df) <- NULL

ggplot(data = filter(mm_MT2_MERVL_df, grepl("MERVL", element_name)), aes(width)) +
  geom_histogram(binwidth = 10, aes(fill = strain)) + 
#     scale_x_continuous(limits = c(4000, 6000)) +
  facet_grid(strain ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
