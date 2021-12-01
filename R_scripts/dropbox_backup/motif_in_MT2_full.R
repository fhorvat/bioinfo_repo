library(dplyr)
library(readr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)

setwd("/common/WORK/fhorvat/Projekti/Svoboda/Other/MT_elements/review/motif_discovery/MT2_full/homer")

# read data
MT2_full <-
  read_csv("../MT2_full_FilteredforCoverage.csv") %>%
  dplyr::select(-repName) %>% 
  dplyr::arrange(FPKM) %>%
  dplyr::slice(1:100) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

# get sequence, write as .FASTA
MT2_full_seq <- as.character(getSeq(x = Mmusculus, MT2_full))

# write FASTA
write.fasta(as.list(MT2_full_seq),
            nbchar = 80,
            names = paste0(MT2_full$fullName, "|", round(MT2_full$FPKM, 3)),
            as.string = TRUE,
            file.out = "MT2_full_FilteredforCoverage_bottom100.fasta",
            open = "w")